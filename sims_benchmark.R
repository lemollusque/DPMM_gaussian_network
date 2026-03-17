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
iter <- 30
dual <- FALSE
param_grid <- expand.grid(
  N = c(100, 200, 500, 1000),
  n = 4:10,
  bge.par = 1
)
sim_grid <- expand.grid(
  j = seq_len(nrow(param_grid)),
  i = seq_len(iter)
)
bge.mus <- c(0.01, 0.1, 0.5, 2, 5)

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
    
    g <- er_dag(n, d=0.2)
    g <- sf_out(g)
    truegraph <- randomize_graph(g)
    model <- cov(truegraph)
    X <- rmvnorm(N, sigma =  model$S)
    data <- standardize(X)
    
    iter_results <- data.frame()
    
    for(k in 1:length(bge.mus)) {
      bge.par <- bge.mus[k]
      bge.searchspace <- set.searchspace(
        data, dual, "bge", bgepar = list(am = bge.par)
      )
      
      bge.fit <- bge.partition.mcmc(bge.searchspace, order = FALSE)
      iter_results <- compare_results(
        bge.fit, c(bge.par, "BGe, partition"), iter_results, truegraph
      )
      
      bge.fit <- bge.partition.mcmc(bge.searchspace, order = TRUE)
      iter_results <- compare_results(
        bge.fit, c(bge.par, "BGe, order"), iter_results, truegraph
      )
    }
    
    iter_results$N <- N
    iter_results$n <- n
    
    iter_results
  }
})
future::plan(sequential)

colnames(results) <- c(
  "ESHD", "eTP", "eFP", "TPR", "FPR_P",
  "time", "parameter", "method", "graph", "N", "n"
)

saveRDS(results, "Results/Sims_benchmark.rds")
#results <- as.data.frame(readRDS("Results/Sims_benchmark.rds"))


# plots resutts
keep_methods <- c("BGe, partition", "BGe, order")

results_small <- results %>%
  filter(graph == "pattern", method %in% keep_methods) 

results_small <- results_small %>%
  mutate(ESHD = as.numeric(ESHD))


best_params <- results_small %>%
  group_by(N, n, method, parameter) %>%
  summarise(mean_ESHD = mean(ESHD, na.rm = TRUE), .groups = "drop") %>%
  group_by(N, n, method) %>%
  slice_min(mean_ESHD, n = 1, with_ties = FALSE)


results_small <- results_small %>%
  inner_join(best_params, by = c("N", "n", "method", "parameter"))


medians <- results_small %>%
  group_by(method, N, n) %>%
  summarise(median_ESHD = median(ESHD), .groups = "drop")

ggplot(results_small, aes(x = method, y = ESHD, color = method)) +
  geom_boxplot(aes(group = method), width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.12, alpha = 0.5) +
  geom_text(
    data = medians,
    aes(x = method, y = median_ESHD, label = round(median_ESHD,2)),
    color = "black",
    vjust = -0.7,
    size = 3
  ) +
  facet_grid(n ~ N, scales = "free_y") +
  theme_bw()
