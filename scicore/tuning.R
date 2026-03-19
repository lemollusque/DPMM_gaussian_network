packages <- c(
  "BiDAG",
  "dirichletprocess",
  "dplyr",
  "ggplot2",
  "foreach",
  "doFuture",
  "future",
  "aricode",
  "matrixStats",
  "progressr",
  "doRNG",
  "mvtnorm"
)

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

source("dao.R")
source("fns.R")

# --------------------------------------------------
# Settings
# --------------------------------------------------

iter <- 10
n <- 4
shift_sizes <- c(1, 3, 5)

N_samples <- c(100, 200, 300, 400, 500)
dp_fits <- 4
dp_iter <- 500

sample_windows <- list(
  w40_50   = 40:50,
  w90_100  = 90:100,
  w190_200 = 190:200,
  w290_300 = 290:300,
  w390_400 = 390:400,
  w490_500 = 490:500
)

alpha_grid <- list(
  c(2, 4)
)

g0_grid <- list(
  list(mu0 = rep(0, n), kappa0 = n, nu = n, Lambda = diag(n))
)

# Hyperparameter grid
param_grid <- expand.grid(
  alpha_id = seq_along(alpha_grid),
  g0_id = seq_along(g0_grid),
  N = N_samples,
  shift_size = shift_sizes
)

# Simulation grid = hyperparameter settings x replicate
sim_grid <- expand.grid(
  hp_id = seq_len(nrow(param_grid)),
  sim_id = seq_len(iter)
)

# --------------------------------------------------
# Helpers
# --------------------------------------------------

dp_membership_probs_iters <- function(dp, sample_iters) {
  y <- dp$data
  N <- dp$n
  mdObj <- dp$mixingDistribution
  
  clusterParams_sample <- dp$clusterParametersChain[sample_iters]
  weightsChain_sample <- dp$weightsChain[sample_iters]
  pointsPerCluster_sample <- lapply(weightsChain_sample, function(w) w * N)
  
  probs_list <- vector("list", length(sample_iters))
  
  for (l in seq_along(sample_iters)) {
    clusterParams <- clusterParams_sample[[l]]
    pointsPerCluster <- pointsPerCluster_sample[[l]]
    numLabels <- length(pointsPerCluster)
    
    probs <- matrix(0, nrow = N, ncol = numLabels)
    
    for (i in seq_len(N)) {
      rowp <- pointsPerCluster *
        Likelihood(mdObj, y[i, , drop = FALSE], clusterParams)
      
      rowp[is.na(rowp)] <- 0
      if (all(rowp == 0)) {
        rowp <- rep_len(1, length(rowp))
      }
      
      probs[i, ] <- rowp
    }
    
    col_sums <- colSums(probs)
    probs <- probs[, col_sums > 0, drop = FALSE]
    probs_list[[l]] <- probs / rowSums(probs)
  }
  
  names(probs_list) <- as.character(sample_iters)
  probs_list
}

make_iter_window_map <- function(sample_windows) {
  bind_rows(lapply(names(sample_windows), function(w) {
    data.frame(
      window = w,
      sampled_iter = sample_windows[[w]]
    )
  }))
}

score_one_gamma <- function(Gamma_sample, truth) {
  x1_mass <- colSums(Gamma_sample[truth == "X1", , drop = FALSE])
  x2_mass <- colSums(Gamma_sample[truth == "X2", , drop = FALSE])
  
  cluster_owner <- ifelse(x1_mass >= x2_mass, "X1", "X2")
  
  same_mass <- numeric(nrow(Gamma_sample))
  other_mass <- numeric(nrow(Gamma_sample))
  
  for (i in seq_len(nrow(Gamma_sample))) {
    same_mass[i]  <- sum(Gamma_sample[i, cluster_owner == truth[i], drop = FALSE])
    other_mass[i] <- sum(Gamma_sample[i, cluster_owner != truth[i], drop = FALSE])
  }
  
  entropy <- -rowSums(Gamma_sample * log(pmax(Gamma_sample, 1e-12)))
  entropy <- entropy / log(ncol(Gamma_sample))
  hard_labels <- max.col(Gamma_sample, ties.method = "first")
  
  data.frame(
    ari = adjustedRandIndex(hard_labels, truth),
    nmi = aricode::NMI(hard_labels, truth),
    same_mass = mean(same_mass),
    other_mass = mean(other_mass),
    entropy = mean(entropy),
    score = mean(same_mass) - mean(other_mass)
  )
}

summarise_hyperparam <- function(df_sim) {
  data.frame(
    n_sim = nrow(df_sim),
    hp_mean_score = mean(df_sim$mean_score),
    hp_sd_score = sd(df_sim$mean_score),
    hp_mean_ari = mean(df_sim$mean_ari),
    hp_sd_ari = sd(df_sim$mean_ari),
    hp_mean_nmi = mean(df_sim$mean_nmi),
    hp_sd_nmi = sd(df_sim$mean_nmi),
    hp_mean_entropy = mean(df_sim$mean_entropy),
    hp_mean_same_mass = mean(df_sim$mean_same_mass),
    hp_mean_other_mass = mean(df_sim$mean_other_mass),
    hp_best_score = max(df_sim$max_score),
    hp_worst_score = min(df_sim$min_score)
  )
}

# --------------------------------------------------
# One simulation
# --------------------------------------------------
run_one_simulation <- function(sim_id, alpha_prior, g0_prior,
                               n, N, dp_fits, dp_iter,
                               sample_windows, shift_size,
                               job_seed) {
  set.seed(job_seed)
  
  sample_iters <- sort(unique(unlist(sample_windows)))
  iter_window_map <- make_iter_window_map(sample_windows)
  
  g <- er_dag(n)
  g <- sf_out(g)
  truegraph <- randomize_graph(g)
  
  model1 <- corr(truegraph)
  model2 <- corr(truegraph)
  
  X1 <- simulate(model1$B, model1$O, N / 2)
  X2 <- simulate(model2$B, model2$O, N / 2)
  
  v <- rnorm(ncol(X2))
  v <- v / sqrt(sum(v^2))
  shift <- shift_size * v
  X2 <- sweep(X2, 2, shift, "+")
  
  data <- standardize(rbind(X1, X2))
  colnames(data) <- paste0("v", seq_len(ncol(data)))
  truth <- c(rep("X1", N / 2), rep("X2", N / 2))
  
  gamma_scores_all <- vector("list", dp_fits * length(sample_iters))
  gamma_counter <- 1
  
  for (f in seq_len(dp_fits)) {
    dp <- DirichletProcessMvnormal(
      data,
      g0Priors = g0_prior,
      alphaPriors = alpha_prior
    )
    
    dp <- Fit(dp, dp_iter)
    gamma_list <- dp_membership_probs_iters(dp, sample_iters)
    
    if (is.matrix(gamma_list)) {
      gamma_list <- list(gamma_list)
    }
    
    for (l in seq_along(gamma_list)) {
      Gamma_sample <- gamma_list[[l]]
      sc <- score_one_gamma(Gamma_sample, truth)
      
      gamma_scores_all[[gamma_counter]] <- data.frame(
        sim = sim_id,
        dp_fit = f,
        gamma_id = l,
        sampled_iter = sample_iters[l],
        ari = sc$ari,
        nmi = sc$nmi,
        same_mass = sc$same_mass,
        other_mass = sc$other_mass,
        entropy = sc$entropy,
        score = sc$score,
        job_seed = job_seed
      )
      gamma_counter <- gamma_counter + 1
    }
  }
  
  gamma_scores_df <- bind_rows(gamma_scores_all) %>%
    left_join(iter_window_map, by = "sampled_iter")
  
  sim_summary_by_window <- gamma_scores_df %>%
    group_by(sim, window, job_seed) %>%
    summarise(
      n_gamma = n(),
      mean_score = mean(score),
      sd_score = sd(score),
      min_score = min(score),
      median_score = median(score),
      max_score = max(score),
      mean_ari = mean(ari),
      mean_nmi = mean(nmi),
      mean_same_mass = mean(same_mass),
      mean_other_mass = mean(other_mass),
      mean_entropy = mean(entropy),
      .groups = "drop"
    )
  
  list(
    sim_summary_by_window = sim_summary_by_window,
    gamma_scores = gamma_scores_df
  )
}


# --------------------------------------------------
# Parallel setup
# --------------------------------------------------

n_cores <- max(1, availableCores() - 1)
plan(multisession, workers = n_cores)
registerDoFuture()

handlers(global = TRUE)
handlers("progress")

# --------------------------------------------------
# Parallel run
# --------------------------------------------------

with_progress({
  p <- progressor(steps = nrow(sim_grid))
  
  foreach(
    k = seq_len(nrow(sim_grid)),
    .packages = c(
      "BiDAG", "matrixStats", "dirichletprocess", "dplyr",
      "aricode", "mvtnorm"
    )
  ) %dopar% {
    
    source("dao.R")
    source("fns.R")
    insertSource("fns.R", package = "BiDAG")
    
    hp_id <- sim_grid$hp_id[k]
    sim_id <- sim_grid$sim_id[k]
    
    alpha_id <- param_grid$alpha_id[hp_id]
    g0_id <- param_grid$g0_id[hp_id]
    N <- param_grid$N[hp_id]
    shift_size <- param_grid$shift_size[hp_id]
    
    alpha_prior <- alpha_grid[[alpha_id]]
    g0_prior <- g0_grid[[g0_id]]
    
    # ------------------------------------
    # unique job seed
    # ------------------------------------
    
    job_seed <- as.integer(
      (as.numeric(Sys.time()) * 1000 + Sys.getpid() + k) %% .Machine$integer.max
    )
    
    if (job_seed <= 0) job_seed <- k
    
    out <- run_one_simulation(
      sim_id = sim_id,
      alpha_prior = alpha_prior,
      g0_prior = g0_prior,
      n = n,
      N = N,
      dp_fits = dp_fits,
      dp_iter = dp_iter,
      sample_windows = sample_windows,
      shift_size = shift_size,
      job_seed = job_seed
    )
    
    sim_summary <- out$sim_summary_by_window %>%
      mutate(
        hp_id = hp_id,
        alpha_a = alpha_prior[1],
        alpha_b = alpha_prior[2],
        g0_id = g0_id,
        kappa0 = g0_prior$kappa0,
        nu = g0_prior$nu,
        lambda_scale = g0_prior$Lambda[1, 1],
        N = N,
        n = n,
        dp_iter = dp_iter,
        shift_size = shift_size,
        job_seed = job_seed
      )
    
    gamma_scores <- out$gamma_scores %>%
      mutate(
        hp_id = hp_id,
        alpha_a = alpha_prior[1],
        alpha_b = alpha_prior[2],
        g0_id = g0_id,
        kappa0 = g0_prior$kappa0,
        nu = g0_prior$nu,
        lambda_scale = g0_prior$Lambda[1, 1],
        N = N,
        n = n,
        dp_iter = dp_iter,
        shift_size = shift_size,
        job_seed = job_seed
      )
    
    result <- bind_rows(
      sim_summary %>% mutate(result_type = "sim_summary"),
      gamma_scores %>% mutate(result_type = "gamma_scores")
    )
    
    # ------------------------------------
    # unique filename
    # ------------------------------------
    
    uid <- paste0(
      format(Sys.time(), "%Y%m%d_%H%M%S"),
      "_pid", Sys.getpid(),
      "_k", k,
      "_seed", job_seed
    )
    
    file_name <- paste0("Tuning/run_", uid, ".rds")
    
    saveRDS(result, file_name)
    
    p()
  }
})
