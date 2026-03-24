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
source("Fourier_fns.R")
source("dao.R")
source("fns.R")
insertSource("fns.R", package = "BiDAG")

init.seed <- 100
iter <- 100
dual <- TRUE

# dirichlet params
dp_iter <- 200
dp_fits <- 2
burnin <- 180
L <- 10

param_grid <- expand.grid(
  N = c(600),
  n = 10,
  d = c(0, 0.5, 1),
  bge.par = 0.01
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
    source("Fourier_fns.R")
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
    
    myDAG <- pcalg::randomDAG(n, prob = 0.2, lB = 1, uB = 2) 
    trueDAG <- as(myDAG, "matrix")
    truegraph <- 1*(trueDAG != 0)
    data <- Fou_nldata(truegraph, N, lambda = d, noise.sd = 1, standardize = T)
    
    if (is.null(colnames(data))) {
      colnames(data) <- paste0("v", seq_len(ncol(data)))
    }
    
    
    # dp settings
    dp_usrpar <- list(
      pctesttype = "bge",
      am = bge.par,
      dp_params = list(dp_iter = dp_iter,
                       dp_fits = dp_fits,
                       burnin = burnin,
                       L = L)
    )
    
    # search spaces
    DP.searchspace <- set.searchspace.fullspace(data, dual, "DP", usrpar = dp_usrpar)
    bge.searchspace <- set.searchspace.fullspace(data, dual, "bge", bge.par)
    
    iter_results <- data.frame()
    
    
    bge.fit <- bge.partition.mcmc(bge.searchspace, order = FALSE)
    iter_results <- compare_results(
      bge.fit, c(bge.par, "BGe, partition"), iter_results, truegraph
    )
    
    bge.fit <- bge.partition.mcmc(bge.searchspace, order = TRUE)
    iter_results <- compare_results(
      bge.fit, c(bge.par, "BGe, order"), iter_results, truegraph
    )
    
    dp.fit <- DP.partition.mcmc(DP.searchspace, order = FALSE)
    iter_results <- compare_results(
      dp.fit, c(bge.par, "DP, partition"), iter_results,  truegraph
    )
    
    dp.fit <- DP.partition.mcmc(DP.searchspace, order = TRUE)
    iter_results <- compare_results(
      dp.fit, c(bge.par, "DP, order"), iter_results,  truegraph
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