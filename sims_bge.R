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
  N = c(1000),
  n = c(10),
  d = c(0, 0.5, 1),
  bge.par = 1
)
sim_grid <- expand.grid(
  j = seq_len(nrow(param_grid)),
  i = seq_len(iter)
)

make_job_seed <- function(init.seed, k) {
  init.seed + k
}

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
  ) %dopar% {
    
    source("comparison_algs.R")
    source("dualPC.R")
    source("dao.R")
    source("fns.R")
    source("Fourier_fns.R")
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
    
    
    
    # deterministic per-job seed
    job_seed <- make_job_seed(init.seed, i +  100)
    set.seed(job_seed)
    
    
    g <- er_dag(n)
    g <- sf_out(g)
    g <- randomize_graph(g)
    truegraph = t(g)
    data <- Fou_nldata(truegraph, N, lambda = d, noise.sd = 1, standardize = T)
    
    if (is.null(colnames(data))) {
      colnames(data) <- paste0("v", seq_len(ncol(data)))
    }
    
    bge.searchspace = set.searchspace(data, dual, "bge", bge.par)
    
    
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

saveRDS(results, "Results/Sims_bimodal_bge.rds")
#results <- as.data.frame(readRDS("Results/Sims_bimodal_bge.rds"))

# plots resutts
keep_methods <- c("BGe, partition", "BGe, order")

results_small <- results %>%
  filter(graph == "pattern", 
         method %in% keep_methods,
         n == 10) 

results_small <- results_small %>%
  mutate(ESHD = as.numeric(ESHD))

medians <- results_small %>%
  group_by(method, N, d) %>%
  summarise(medians_ESHD = median(ESHD), .groups = "drop")

ggplot(results_small, aes(x = method, y = ESHD, color = method)) +
  geom_boxplot(aes(group = method), width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.12, alpha = 0.5) +
  geom_text(
    data = medians,
    aes(x = method, y = medians_ESHD, label = round(medians_ESHD,2)),
    color = "black",
    vjust = -0.7,
    size = 3
  ) +
  facet_grid(d ~ N, scales = "free_y") +
  theme_bw()
