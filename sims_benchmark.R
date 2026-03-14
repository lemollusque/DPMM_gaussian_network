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

source("comparison_algs.R")
source("dualPC.R")
source("dao.R")
source("fns.R")
insertSource("fns.R", package = "BiDAG")

init.seed <- 100
iter <- 100
dual <- FALSE
param_grid <- expand.grid(
  N = c(100,200),
  n = 4,
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
    bge.par <- param_grid$bge.par[j]
    
    g <- er_dag(n, ad = 3)
    g <- sf_out(g)
    truegraph <- randomize_graph(g)
    
    model <- corr(g)
    
    df <- 3
    t_err <- function(n, var) {
      rt(n, df = df) * sqrt(var * (df - 2) / df)
    }
    
    X <- simulate(model$B, model$O, N)
    data <- standardize(X)
    
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
    
    iter_results
  }
})

colnames(results) <- c(
  "ESHD", "eTP", "eFP", "TPR", "FPR_P",
  "time", "parameter", "method", "graph", "N", "n"
)

saveRDS(results, "Results/Sims_benchmark.rds")


# plots resutts
keep_methods <- c("BGe, partition", "BGe, order")

results_small <- results %>%
  filter(graph == "pattern", method %in% keep_methods) 

results_small <- results_small %>%
  mutate(ESHD = as.numeric(ESHD))

ggplot(results_small, aes(x = method, y = ESHD, color = method)) +
  geom_boxplot(aes(group = method), width = 0.6, outlier.shape = NA, linewidth = 0.6) +
  geom_jitter(aes(group = method), width = 0.15, alpha = 0.7, size = 1.2) +
  coord_cartesian(ylim = c(0, 12)) +
  labs(x = NULL, y = "E=SHD") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major.x = element_blank()
  )
