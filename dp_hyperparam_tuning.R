library(BiDAG)
library(matrixStats)
library(dirichletprocess)
library(dplyr)
library(ggplot2)
library(foreach)
library(doFuture)
library(future)
library(parallelly)
library(aricode)
library(mclust)
library(openxlsx)
library(progressr)
library(mvtnorm)
library(doRNG)

source("dao.R")
source("dualPC.R")
source("Fourier_fns.R")
source("fns.R")

# --------------------------------------------------
# Settings
# --------------------------------------------------

init.seed <- 100
iter <- 100
dual <- TRUE

n <- 10
shift_sizes <- c(10)
N_samples <- c(100)

dp_fits <- 1
dp_iter <- 100

alpha_grid <- list(
  c(2,4)
)

g0_grid <- list(
  list(mu0 = rep(0, n), kappa0 = 1, nu = n, Lambda = diag(n)/(n))
)

# --------------------------------------------------
# Hyperparameter grid
# --------------------------------------------------

param_grid <- expand.grid(
  alpha_id = seq_along(alpha_grid),
  g0_id = seq_along(g0_grid),
  N = N_samples,
  shift_size = shift_sizes
)

sim_grid <- expand.grid(
  j = seq_len(nrow(param_grid)),
  i = seq_len(iter)
)

# --------------------------------------------------
# Helpers
# --------------------------------------------------

dp_nclusters_iters <- function(dp, step = 25) {
  weights_sample <- dp$weightsChain
  
  n_clusters <- sapply(weights_sample, function(w) {
    sum(w > 0)
  })
  
  iters <- seq_along(n_clusters)
  keep <- iters %% step == 0 | iters == 1
  
  data.frame(
    iter = iters[keep],
    n_clusters = as.integer(n_clusters[keep])
  )
}

make_job_seed <- function(init.seed, k) {
  init.seed + k
}

make_file_name <- function(alpha_id, kappa, lambda_coef, shift_size, N, i) {
  paste0(
    "Tuning/",
    "alpha_id", alpha_id,
    "_g0_id", kappa,"_", lambda_coef,
    "_shift_size", shift_size,
    "_N", N,
    "_rep", sprintf("%03d", i),
    ".rds"
  )
}

# --------------------------------------------------
# Parallel setup
# --------------------------------------------------
n_cores <- max(1, availableCores() - 1)

plan(multisession, workers = n_cores)
registerDoFuture()
registerDoRNG(init.seed)

handlers(global = TRUE)
handlers("txtprogressbar")

# --------------------------------------------------
# Parallel run
# --------------------------------------------------

with_progress({
  p <- progressor(steps = nrow(sim_grid))
  
  foreach(
    k = seq_len(nrow(sim_grid)),
    .packages = c(
      "BiDAG", "matrixStats", "dirichletprocess", "dplyr",
      "aricode", "mclust", "mvtnorm", "doRNG", "pcalg"
    )
  ) %dopar% {
    
    source("dao.R")
    source("dualPC.R")
    source("Fourier_fns.R")
    source("fns.R")
    insertSource("fns.R", package = "BiDAG")
    
    # set params
    j <- sim_grid$j[k]
    i <- sim_grid$i[k]
    
    alpha_id <- param_grid$alpha_id[j]
    g0_id <- param_grid$g0_id[j]
    
    N <- param_grid$N[j]
    shift_size <- param_grid$shift_size[j]
    
    alpha_prior <- alpha_grid[[alpha_id]]
    alpha = alpha_prior[[1]]/alpha_prior[[2]]
    g0_prior <- g0_grid[[g0_id]]
    kappa = g0_prior$kappa0
    lambda_coef = g0_prior$Lambda[1,1]
    
    file_name <- make_file_name(
      alpha_id = alpha,
      kappa = kappa,
      lambda_coef = lambda_coef,
      shift_size = shift_size,
      N = N,
      i = i
    )
    
    if (file.exists(file_name)) {
      p(sprintf("skip %d", k))
      return(NULL)
    }
    
    job_seed <- make_job_seed(init.seed, k)
    set.seed(job_seed)
    
    myDAG <- pcalg::randomDAG(n, prob = 0.6, lB = 1, uB = 2)
    trueDAG <- as(myDAG, "matrix")
    truegraph <- 1 * (trueDAG != 0)
    
    data <- simulate_bimodal(truegraph, n=N, bimodal_sep=shift_size)
    
    if (is.null(colnames(data))) {
      colnames(data) <- paste0("v", seq_len(ncol(data)))
    }
    
    iter_results <- data.frame()
    
    # full DP
    for (f in seq_len(dp_fits)) {
      dp <- DirichletProcessMvnormal(
        data,
        g0Priors = g0_prior,
        alphaPriors = alpha_prior,
        numInitialClusters = nrow(data)
      )
      dp <- Fit(dp, dp_iter, progressBar = FALSE)
      
      fit_results <- dp_nclusters_iters(dp)
      fit_results$dp_fit <- f
      
      iter_results <- bind_rows(iter_results, fit_results)
    }
    iter_results$alpha <- alpha
    iter_results$kappa <- kappa
    iter_results$lambda_coef <- lambda_coef
    iter_results$shift_size <- shift_size
    iter_results$N <- N
    iter_results$rep <- i
    
    saveRDS(iter_results, file_name)
    
    p(sprintf("done %d", k))
  }
})

plan(sequential)

# --------------------------------------------------
# load combined runs
# --------------------------------------------------

files <- list.files("Tuning", pattern = "\\.rds$", full.names = TRUE)
results <- bind_rows(lapply(files, readRDS))

# --------------------------------------------------
# Summaries
# --------------------------------------------------

summary_results <- results %>%
  group_by(alpha, kappa, lambda_coef, N, shift_size, iter) %>%
  summarise(
    mean_clusters = mean(n_clusters, na.rm = TRUE),
    max_clusters = max(n_clusters, na.rm = TRUE),
    sd_clusters = sd(n_clusters, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(kappa_lambda = paste0("k=", kappa, ", T=", lambda_coef))

label_df <- summary_results %>%
  filter(iter %% 50 == 0) %>%
  mutate(label = round(mean_clusters, 2))

# --------------------------------------------------
# Plot
# --------------------------------------------------

ggplot() +
  geom_point(
    data = results,
    aes(x = iter, y = n_clusters),
    alpha = 0.4
  ) +
  geom_line(
    data = summary_results,
    aes(x = iter, y = mean_clusters, group = 1),
    linewidth = 1
  ) +
  geom_text(
    data = label_df,
    aes(x = iter, y = mean_clusters, label = label),
    vjust = -0.5,
    size = 3
  )+
  facet_grid(kappa_lambda ~ alpha) +
  theme_bw()