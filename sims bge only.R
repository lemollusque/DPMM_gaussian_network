library(BiDAG)
library(matrixStats)
library(foreach)
library(doFuture)
library(future)
library(parallelly)
library(aricode)
library(mclust)
library(openxlsx)
library(progressr)
library(doRNG)
library(tidyverse)
library(mvtnorm)

source("Fourier_fns.R")
source("BayesStanFns.R")
source("sampling_fns.R")
source("comparison_algs.R")
source("dualPC.R")
source("dao.R")
insertSource("GPscore.R", package = "BiDAG")

init.seed <- 100
iter <- 30  # number of simulations
lambdas <- c(0)  # non-linearity; zero is linear
dual <- T    # use dualPC
n <- 10      # number of nodes
N <- 100     # number of samples
results <- data.frame()

# Parameters for ROC curves
bge.mus <- c(0.01, 0.1, 0.5, 2, 5)

# Grid over all independent simulation jobs
sim_grid <- expand.grid(
  lambda = lambdas,
  bge.par = bge.mus,
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
    .packages = c("BiDAG", "matrixStats", "dplyr")
  ) %dorng% {
    
    # Re-source on each worker
    source("Fourier_fns.R")
    source("BayesStanFns.R")
    source("sampling_fns.R")
    source("comparison_algs.R")
    source("dualPC.R")
    source("dao.R")
    insertSource("GPscore.R", package = "BiDAG")
    
    p()
    lambda <- sim_grid$lambda[k]
    bge.par <- sim_grid$bge.par[k]
    i <- sim_grid$i[k]
    # optional: explicit seed per task
    set.seed(init.seed + i)
    # Generate DAG & data
    myDAG <- pcalg::randomDAG(n, prob = 0.2, lB = 1, uB = 2) 
    trueDAG <- as(myDAG, "matrix")
    truegraph <- 1*(trueDAG != 0)
    
    model <- corr(truegraph)
    
    X <- simulate(model$B, model$O, N)
    data <- standardize(X)
    
    
    iter_results <- data.frame()
    
    bge.searchspace = set.searchspace(data, dual, "bge", bge.par)
    
    # Bge score, partition
    bge.fit <- bge.partition.mcmc(bge.searchspace, order = F)
    iter_results <- compare_results(bge.fit, c(bge.par, "BGe, partition", lambda), iter_results, truegraph)
    
    # Bge score, order
    bge.fit <- bge.partition.mcmc(bge.searchspace, order = T)
    iter_results <- compare_results(bge.fit, c(bge.par, "BGe, order", lambda), iter_results, truegraph)
    # add identifiers
    iter_results$N <- N
    iter_results$n <- n
    iter_results$iter <- i
    
    iter_results
  }
})
future::plan(sequential)

colnames(results) <- c("ESHD", "eTP", "eFP", "TPR", "FPR_P", 
                       "time", "parameter", "method", "lambda", "graph",
                       "N", "n", "iter"
)
saveRDS(results, "Results/Sims_Results_bge.rds")

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
