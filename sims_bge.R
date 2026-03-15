library(BiDAG)
library(matrixStats)
library(dirichletprocess)
library(dplyr)
library(ggplot2)
library(foreach)
library(doFuture)
library(future)
library(parallelly)
library(mclust)
library(progressr)
library(doRNG)
library(mvtnorm)

source("comparison_algs.R")
source("dualPC.R")
source("dao.R")
source("fns.R")
insertSource("fns.R", package = "BiDAG")

init.seed <- 100
iter <- 100
dual <- FALSE
param_grid <- expand.grid(
  N = c(100, 200, 500, 1000),
  n = 4:10,
  d = 0:10,
  bge.par = 1
)
sim_grid <- expand.grid(
  j = seq_len(nrow(param_grid)),
  i = seq_len(iter)
)

n_cores <- max(1, availableCores() - 1)

plan(multisession, workers = n_cores)
registerDoFuture()
registerDoRNG(init.seed)

handlers(global = TRUE)
handlers("txtprogressbar")

results <- with_progress({
  p <- progressor(steps = nrow(sim_grid))
  
  foreach(
    k = seq_len(nrow(sim_grid)),
    .combine = rbind,
    .packages = c("BiDAG", "matrixStats", "dirichletprocess", "dplyr", "mclust")
  ) %dorng% {
    
    source("comparison_algs.R")
    source("dualPC.R")
    source("dao.R")
    source("fns.R")
    insertSource("fns.R", package = "BiDAG")
    
    # show progress
    p()
    
    # set params
    j <- sim_grid$j[k]
    i <- sim_grid$i[k]
  
    N <- param_grid$N[j]
    n <- param_grid$n[j]
    d <- param_grid$d[j]
    bge.par <- param_grid$bge.par[j]
    
    g <- er_dag(n)
    g <- sf_out(g)
    truegraph <- randomize_graph(g)
    
    model1 <- cov(g)
    model2 <- cov(g)
    
    
    X1 <- rmvt(N/2, sigma =  model1$S, df = 3)
    X2 <- rmvt(N/2, sigma =  model2$S, df = 3)
    
    v <- rnorm(ncol(X2))
    v <- v / sqrt(sum(v^2))   
    shift <- d * v
    X2 <- sweep(X2, 2, shift, "+")
    
    data <- standardize(rbind(X1, X2))
    
    bge.searchspace <- set.searchspace(
      data, dual, "bge", bgepar = list(am = bge.par)
    )
    
    iter_results <- data.frame()
    
    bge.fit <- bge.partition.mcmc(bge.searchspace, order = FALSE)
    iter_results <- compare_results(
      bge.fit, c(bge.par, "BGe, partition"), iter_results, truegraph
    )
    
    bge.fit <- bge.partition.mcmc(bge.searchspace, order = TRUE)
    iter_results <- compare_results(
      bge.fit, c(bge.par, "BGe, order"), iter_results, truegraph
    )
    
    iter_results$N <- N
    iter_results$n <- n
    iter_results$d <- d
    
    iter_results
  }
})
future::plan(sequential)

colnames(results) <- c(
  "ESHD", "eTP", "eFP", "TPR", "FPR_P",
  "time", "parameter", "method", "graph", "N", "n", "d"
)

saveRDS(results, "Results/Sims_bimodal_bge_student.rds")


# plots resutts
keep_methods <- c("BGe, partition", "BGe, order")

results_small <- results %>%
  filter(graph == "pattern", 
         method %in% keep_methods,
         n == 10) 

results_small <- results_small %>%
  mutate(ESHD = as.numeric(ESHD))

means <- results_small %>%
  group_by(method, N, d) %>%
  summarise(mean_ESHD = mean(ESHD), .groups = "drop")

ggplot(results_small, aes(x = method, y = ESHD, color = method)) +
  geom_boxplot(aes(group = method), width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.12, alpha = 0.5) +
  geom_text(
    data = means,
    aes(x = method, y = mean_ESHD, label = round(mean_ESHD,2)),
    color = "black",
    vjust = -0.7,
    size = 3
  ) +
  facet_grid(d ~ N, scales = "free_y") +
  theme_bw()
