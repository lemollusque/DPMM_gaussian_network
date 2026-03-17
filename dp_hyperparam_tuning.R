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
library(doRNG)
library(mvtnorm)

source("dao.R")
source("fns.R")

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

summarise_scores <- function(df_scores) {
  data.frame(
    n_gamma = nrow(df_scores),
    mean_score = mean(df_scores$score),
    sd_score = sd(df_scores$score),
    min_score = min(df_scores$score),
    median_score = median(df_scores$score),
    max_score = max(df_scores$score),
    mean_ari = mean(df_scores$ari),
    mean_nmi = mean(df_scores$nmi),
    mean_same_mass = mean(df_scores$same_mass),
    mean_other_mass = mean(df_scores$other_mass),
    mean_entropy = mean(df_scores$entropy)
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

run_one_simulation <- function(sim_id, alpha_prior, g0_prior,
                               init.seed, n, N, dp_fits, dp_iter,
                               burnin, L, shift_size) {
  set.seed(init.seed + sim_id)
  
  g <- er_dag(n, d=0.2)
  g <- sf_out(g)
  truegraph <- randomize_graph(g)
  
  model1 <- cov(truegraph)
  model2 <- cov(truegraph)
  
  X1 <- rmvt(N / 2, sigma = model1$S, df = 3)
  X2 <- rmvt(N / 2, sigma = model2$S, df = 3)
  
  v <- rnorm(ncol(X2))
  v <- v / sqrt(sum(v^2))
  shift <- shift_size * v
  X2 <- sweep(X2, 2, shift, "+")
  
  data <- standardize(rbind(X1, X2))
  colnames(data) <- paste0("v", seq_len(ncol(data)))
  truth <- c(rep("X1", N / 2), rep("X2", N / 2))
  
  gamma_scores_all <- vector("list", length = dp_fits * L)
  gamma_counter <- 1
  
  for (f in seq_len(dp_fits)) {
    dp <- DirichletProcessMvnormal(
      data,
      g0Priors = g0_prior,
      alphaPriors = alpha_prior
    )
    
    dp <- Fit(dp, dp_iter)
    gamma_list <- dp_membership_probs(dp, burnin, L)
    
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
        ari = sc$ari,
        nmi = sc$nmi,
        same_mass = sc$same_mass,
        other_mass = sc$other_mass,
        entropy = sc$entropy,
        score = sc$score
      )
      gamma_counter <- gamma_counter + 1
    }
  }
  
  gamma_scores_df <- bind_rows(gamma_scores_all)
  sim_summary <- summarise_scores(gamma_scores_df)
  sim_summary$sim <- sim_id
  
  list(
    sim_summary = sim_summary,
    gamma_scores = gamma_scores_df
  )
}

# ---------------------------
# Settings
# ---------------------------
init.seed <- 100
iter <- 7
n <- 10
shift_size <- 4   

N_samples <- c(2000)
dp_fits <- 1
dp_iter_list <- c(100)

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
  dp_iter = dp_iter_list
)

# Simulation grid = hyperparameter settings x replicate
sim_grid <- expand.grid(
  hp_id = seq_len(nrow(param_grid)),
  sim_id = seq_len(iter)
)

# ---------------------------
# Parallel setup
# ---------------------------
n_cores <- max(1, availableCores() - 1)
plan(multisession, workers = n_cores)
registerDoFuture()
registerDoRNG(init.seed)

handlers(global = TRUE)
handlers("progress")
# ---------------------------
# Parallel run
# ---------------------------
all_runs <- with_progress({
  p <- progressor(steps = nrow(sim_grid))
  
  foreach(
    k = seq_len(nrow(sim_grid)),
    .combine = dplyr::bind_rows,
    .packages = c(
      "BiDAG", "matrixStats", "dirichletprocess", "dplyr",
      "aricode", "mclust", "mvtnorm"
    )
  ) %dorng% {
    
    source("dao.R")
    source("fns.R")
    insertSource("fns.R", package = "BiDAG")
    
    hp_id <- sim_grid$hp_id[k]
    sim_id <- sim_grid$sim_id[k]
    
    alpha_id <- param_grid$alpha_id[hp_id]
    g0_id <- param_grid$g0_id[hp_id]
    N <- param_grid$N[hp_id]
    dp_iter <- param_grid$dp_iter[hp_id]
    
    alpha_prior <- alpha_grid[[alpha_id]]
    g0_prior <- g0_grid[[g0_id]]
    
    L <- 10
    burnin <- dp_iter - L
    
    out <- run_one_simulation(
      sim_id = sim_id,
      alpha_prior = alpha_prior,
      g0_prior = g0_prior,
      init.seed = init.seed + 10000 * hp_id,  # helps separate settings
      n = n,
      N = N,
      dp_fits = dp_fits,
      dp_iter = dp_iter,
      burnin = burnin,
      L = L,
      shift_size = shift_size
    )
    
    sim_summary <- out$sim_summary %>%
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
        dp_iter = dp_iter
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
        dp_iter = dp_iter
      )
    
    # show progress
    p()
    
    data.frame(
      result_type = c(rep("sim_summary", nrow(sim_summary)),
                      rep("gamma_scores", nrow(gamma_scores))),
      bind_rows(sim_summary, gamma_scores)
    )
  }
})


plan(sequential)

# ---------------------------
# Split outputs back apart
# ---------------------------
all_sim_results_df <- all_runs %>%
  filter(result_type == "gamma_scores") %>%
  select(-result_type)

sim_level_summary_df <- all_runs %>%
  filter(result_type == "sim_summary") %>%
  select(-result_type)

# Hyperparameter-level summary
results <- sim_level_summary_df %>%
  group_by(
    hp_id, N, n, dp_iter, alpha_a, alpha_b,
    g0_id, kappa0, nu, lambda_scale
  ) %>%
  group_modify(~ summarise_hyperparam(.x)) %>%
  ungroup()

# ---------------------------
# Save
# ---------------------------
saveRDS(results, "Results/hyper_param_results.rds")
saveRDS(sim_level_summary_df, "Results/hyper_param_sim_level_summaries.rds")
saveRDS(all_sim_results_df, "Results/hyper_param_sim_level_results.rds")

write.xlsx(results, "Results/hyper_param_results.xlsx")
write.xlsx(sim_level_summary_df, "Results/hyper_param_sim_level_summaries.xlsx")
write.xlsx(all_sim_results_df, "Results/hyper_param_sim_level_results.xlsx")

# ---------------------------
# Plot
# ---------------------------
ggplot(results, aes(x = N, y = hp_mean_ari)) +
  geom_point() +
  geom_line(aes(group = g0_id)) +
  facet_wrap(~ alpha_a, scales = "free_y") +
  coord_flip()
