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

# dirichlet params
alpha_prior <- c(2, 4)
initial_clusters <- 10
dp_iter <- 100
dp_fits <- 1
burnin <- 90
L <- 10

param_grid <- expand.grid(
  N = c(1000),
  n = 4,
  d = c(0, 1, 2, 5, 10),
  bge.par = 1
)

sim_grid <- expand.grid(
  j = seq_len(nrow(param_grid)),
  i = seq_len(iter)
)

make_job_seed <- function(init.seed, k) {
  init.seed + k
}

make_file_name <- function(N, n, d, bge.par, i) {
  paste0(
    "Sims/",
    "N", N,
    "_n", n,
    "_d", d,
    "_bge", bge.par,
    "_rep", sprintf("%03d", i),
    ".rds"
  )
}

n_cores <- max(1, availableCores() - 1)

plan(multisession, workers = n_cores)
registerDoFuture()
registerDoRNG(init.seed)

handlers(global = TRUE)
handlers("txtprogressbar")

with_progress({
  p <- progressor(steps = nrow(sim_grid))
  
  foreach(
    k = seq_len(nrow(sim_grid)),
    .packages = c("BiDAG", "matrixStats", "dirichletprocess", "dplyr", "mclust", "mvtnorm")
  ) %dopar% {
    
    source("comparison_algs.R")
    source("dualPC.R")
    source("dao.R")
    source("fns.R")
    insertSource("fns.R", package = "BiDAG")
    
    # set params
    j <- sim_grid$j[k]
    i <- sim_grid$i[k]
    
    N <- param_grid$N[j]
    n <- param_grid$n[j]
    d <- param_grid$d[j]
    bge.par <- param_grid$bge.par[j]
    
    file_name <- make_file_name(
      N = N,
      n = n,
      d = d,
      bge.par = bge.par,
      i = i
    )
    
    # resume logic: skip completed jobs
    if (file.exists(file_name)) {
      p(sprintf("skip %d", k))
      return(NULL)
    }
    
    # deterministic per-job seed
    job_seed <- make_job_seed(init.seed, i)
    set.seed(job_seed)
    
    
    g <- er_dag(n)
    g <- sf_out(g)
    truegraph <- randomize_graph(g)
    
    model <- corr(truegraph)
    X <- simulate_bimodal(model$B, model$O, n=N, bimodal_sep=d)
    data <- standardize(X)
    if (is.null(colnames(data))) {
      colnames(data) <- paste0("v", seq_len(ncol(data)))
    }
    
    # prepare dirichlet gamma list
    Gamma_list <- list()
    for (f in seq_len(dp_fits)) {
      dp <- DirichletProcessMvnormal(data, 
                                     alphaPriors = alpha_prior,
                                     numInitialClusters = initial_clusters)
      dp <- Fit(dp, dp_iter, progressBar = FALSE)
      
      Gamma_sample <- dp_membership_probs(dp, burnin, L)
      Gamma_list <- add_membershipp(
        Gamma_list,
        Gamma_sample,
        child = colnames(data)[1],
        parents = colnames(data)[2:ncol(data)],
        active = TRUE
      )
    }
    
    dp_usrpar <- list(
      pctesttype = "bge",
      am = bge.par,
      membershipp_list = Gamma_list,
      edgepf = 1
    )
    
    # search spaces
    DP.searchspace <- set.searchspace(data, dual, "DP", usrpar = dp_usrpar)
    bge.searchspace <- set.searchspace(data, dual, "bge", bge.par)
    
    iter_results <- data.frame()
    
    bge.fit <- bge.partition.mcmc(bge.searchspace, order = FALSE)
    iter_results <- compare_results(
      bge.fit, c(bge.par, "BGe, partition"), iter_results,  t(truegraph)
    )
    
    bge.fit <- bge.partition.mcmc(bge.searchspace, order = TRUE)
    iter_results <- compare_results(
      bge.fit, c(bge.par, "BGe, order"), iter_results,  t(truegraph)
    )
    
    dp.fit <- DP.partition.mcmc(DP.searchspace, order = FALSE)
    iter_results <- compare_results(
      dp.fit, c(bge.par, "DP, partition"), iter_results,  t(truegraph)
    )
    
    dp.fit <- DP.partition.mcmc(DP.searchspace, order = TRUE)
    iter_results <- compare_results(
      dp.fit, c(bge.par, "DP, order"), iter_results,  t(truegraph)
    )
    
    iter_results$N <- N
    iter_results$n <- n
    iter_results$d <- d
    iter_results$rep <- i
    iter_results$bge.par <- bge.par
    iter_results$job_id <- k
    iter_results$job_seed <- job_seed
    
    colnames(iter_results) <- c(
      "ESHD", "eTP", "eFP", "TPR", "FPR_P",
      "time", "parameter", "method", "graph",
      "N", "n", "d", "rep", "bge.par", "job_id", "job_seed"
    )
    
    saveRDS(iter_results, file_name)
    
    p(sprintf("done %d", k))
  }
})

future::plan(sequential)

files <- list.files("Sims", pattern = "\\.rds$", full.names = TRUE)
results <- bind_rows(lapply(files, readRDS))

# plots resutts
keep_methods <- c("DP, partition", "DP, order", "BGe, partition", "BGe, order")

results_small <- results %>%
  filter(graph == "pattern", method %in% keep_methods) 

results_small <- results_small %>%
  mutate(ESHD = as.numeric(ESHD))
# add text medians
medians <- results_small %>%
  group_by(method, N, n, d) %>%
  summarise(median_ESHD = median(ESHD), .groups = "drop")

ggplot(results_small, aes(x = method, y = ESHD, color = method)) +
  geom_boxplot(aes(group = method), width = 0.6, outlier.shape = NA, linewidth = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.7, size = 1.2) +
  geom_text(
    data = medians,
    aes(x = method, y = median_ESHD, label = round(median_ESHD,2)),
    color = "black",
    vjust = -0.7,
    size = 3
  ) +
  
  labs(x = NULL, y = "E=SHD") +
  facet_grid( ~ d, scales = "free_y") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major.x = element_blank()
  )