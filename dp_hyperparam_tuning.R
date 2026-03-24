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

init.seed <- 102
iter <- 7
dual <- TRUE

n <- 10
shift_sizes <- c(1)
N_samples <- c(102)

dp_fits <- 1
dp_iter <- 400
initial_clusters <- 10
d <- 1

alpha_grid <- list(
  c(10, 5.3)
)

g0_grid <- list(
  list(mu0 = rep(0, n), kappa0 = 1, nu = n, Lambda = diag(n))
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

dp_nclusters_iters <- function(dp) {
  weights_sample <- dp$weightsChain
  
  n_clusters <- sapply(weights_sample, function(w) {
    sum(w > 0)
  })
  
  data.frame(
    iter = seq_along(n_clusters),
    n_clusters = as.integer(n_clusters)
  )
}

make_job_seed <- function(init.seed, k) {
  init.seed + k
}

make_file_name <- function(alpha_id, g0_id, shift_size, N, i) {
  paste0(
    "Tuning/",
    "alpha_id", alpha_id,
    "_g0_id", g0_id,
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
    g0_prior <- g0_grid[[g0_id]]
    
    file_name <- make_file_name(
      alpha_id = alpha_id,
      g0_id = g0_id,
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
    
    myDAG <- pcalg::randomDAG(n, prob = 0.2, lB = 1, uB = 2)
    trueDAG <- as(myDAG, "matrix")
    truegraph <- 1 * (trueDAG != 0)
    
    data <- Fou_nldata(truegraph, N, lambda = d, noise.sd = 1, standardize = TRUE)
    
    if (is.null(colnames(data))) {
      colnames(data) <- paste0("v", seq_len(ncol(data)))
    }
    
    iter_results <- data.frame()
    
    # full DP
    for (f in seq_len(dp_fits)) {
      dp <- DirichletProcessMvnormal(
        data,
        alphaPriors = alpha_prior,
        numInitialClusters = initial_clusters
      )
      dp <- Fit(dp, dp_iter, progressBar = FALSE)
      
      fit_results <- dp_nclusters_iters(dp)
      fit_results$mode <- "full"
      fit_results$dp_fit <- f
      
      iter_results <- bind_rows(iter_results, fit_results)
    }
    
    # dual DP
    if (dual) {
      alpha <- 0.05
      cor_mat <- cor(data)
      startspace <- dual_pc(cor_mat, nrow(data), alpha = alpha, skeleton = TRUE)
      
      idx <- which(rowSums(startspace) > 0 | colSums(startspace) > 0)
      nodes <- rownames(startspace)[idx]
      
      if (length(nodes) >= 2) {
        for (f in seq_len(dp_fits)) {
          dp <- DirichletProcessMvnormal(
            data[, nodes, drop = FALSE],
            alphaPriors = alpha_prior,
            numInitialClusters = initial_clusters
          )
          dp <- Fit(dp, dp_iter, progressBar = FALSE)
          
          fit_results <- dp_nclusters_iters(dp)
          fit_results$mode <- "dual"
          fit_results$dp_fit <- f
          
          iter_results <- bind_rows(iter_results, fit_results)
        }
      }
    }
    
    iter_results$alpha_id <- alpha_id
    iter_results$g0_id <- g0_id
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
  group_by(mode, N, shift_size, iter) %>%
  summarise(
    mean_clusters = mean(n_clusters, na.rm = TRUE),
    sd_clusters = sd(n_clusters, na.rm = TRUE),
    .groups = "drop"
  )
label_df <- summary_results %>%
  filter(iter %% 50 == 0) %>%
  mutate(label = round(mean_clusters, 2))

# --------------------------------------------------
# Plot
# --------------------------------------------------
ggplot(summary_results, aes(x = iter, y = mean_clusters)) +
  geom_line() +
  geom_point() +
  geom_errorbar(
    aes(
      ymin = mean_clusters - sd_clusters,
      ymax = mean_clusters + sd_clusters
    ),
    width = 3
  ) +
  geom_text(
    data = label_df,
    aes(y = 5, label = label),
    vjust = -0.5,
    size = 3
  ) +
  ylim(c(0, 10)) +
  facet_grid(N ~ mode) +
  theme_bw()
