library(BiDAG)
library(matrixStats)
library(dirichletprocess)
library(dplyr)
library(ggplot2)
library(foreach)
library(doParallel)
library(parallelly)
library(aricode)
library(mclust)

source("dao.R")
source("fns.R")



init.seed <- 100
iter <- 3  # number of simulations per hyperparam
n <- 4    # number of nodes
N <- 100     # number of samples
results <- data.frame()

dp_fits = 2
dp_iter = 100
burnin = 50
L = 10

truth <- c(rep("X1", N/2), rep("X2", N/2))


score_one_gamma <- function(Gamma_sample, truth) {
  # Gamma_sample: N x K matrix
  # truth: length N vector with values "X1"/"X2"
  
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
    hp_mean_entropy = mean(df_sim$mean_entropy),
    hp_mean_same_mass = mean(df_sim$mean_same_mass),
    hp_mean_other_mass = mean(df_sim$mean_other_mass),
    hp_best_score = max(df_sim$max_score),
    hp_worst_score = min(df_sim$min_score)
  )
}

run_one_simulation <- function(sim_id, alpha_prior, g0_prior, 
                               init.seed, n, N, dp_fits, dp_iter, burnin, L) {
  set.seed(init.seed + sim_id)
  
  g <- er_dag(n, ad = 3)
  g <- sf_out(g)
  truegraph <- randomize_graph(g)
  
  model1 <- corr(g)
  model2 <- corr(g)
  
  df <- 3
  t_err <- function(n, var) {
    rt(n, df = df) * sqrt(var * (df - 2) / df)
  }
  
  X1 <- simulate(model1$B, model1$O, N/2, err = t_err)
  X2 <- simulate(model2$B, model2$O, N/2, err = t_err)
  
  shift <- runif(ncol(X2), -2, 2)
  X2 <- sweep(X2, 2, shift, "+")
  
  data <- standardize(rbind(X1, X2))
  colnames(data) <- paste0("v", seq_len(ncol(data)))
  truth <- c(rep("X1", N/2), rep("X2", N/2))
  
  gamma_scores_all <- list()
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
  sim_summary
}

alpha_grid <- list(
  c(1, 1),
  c(2, 4),
  c(2, 2),
  c(4, 2)
)

alpha_mu <- 1          
alpha_w  <- n + alpha_mu + 1      
t <- alpha_mu * (alpha_w - n - 1) / (alpha_mu + 1)

g0_grid <- list(
  list(mu0 = rep(0, n), kappa0 = n,   nu = n,   Lambda = diag(n)),
  list(mu0 = rep(0, n), Lambda = diag(n) / t, kappa0 = alpha_mu, nu = alpha_w)
)

results <- data.frame()
all_sim_results <- list()
counter <- 1

sink("logfile.txt")

for (a in seq_along(alpha_grid)) {
  for (g in seq_along(g0_grid)) {
      
      cat("Hyperparam setting:", counter, "\n")
      cat("alpha =", paste(alpha_grid[[a]], collapse = ","), "\n")
      cat("g0 =", g, "\n")
      
      sim_results <- list()
      
      for (i in seq_len(iter)) {
        cat("  Simulation:", i, "\n")
        
        sim_results[[i]] <- run_one_simulation(
          sim_id = i,
          alpha_prior = alpha_grid[[a]],
          g0_prior = g0_grid[[g]],
          init.seed = init.seed,
          n = n,
          N = N,
          dp_fits = dp_fits,
          dp_iter = dp_iter,
          burnin = burnin,
          L = L
        )
      }
      
      sim_results_df <- bind_rows(sim_results)
      hp_summary <- summarise_hyperparam(sim_results_df)
      
      hp_summary$alpha_a <- alpha_grid[[a]][1]
      hp_summary$alpha_b <- alpha_grid[[a]][2]
      hp_summary$g0_id <- g
      
      results <- bind_rows(results, hp_summary)
      
      all_sim_results[[counter]] <- sim_results_df %>%
        mutate(
          alpha_a = alpha_grid[[a]][1],
          alpha_b = alpha_grid[[a]][2],
          g0_id = g,
        )
      
      counter <- counter + 1
  }
}

sink()

all_sim_results_df <- bind_rows(all_sim_results)

saveRDS(results, "Results/hyper_param_results.rds")
saveRDS(all_sim_results_df, "Results/hyper_param_sim_level_results.rds")

ggplot(results, aes(x = interaction(alpha_a, alpha_b, g0_id),
                    y = hp_mean_ari)) +
  geom_point() +
  coord_flip()
