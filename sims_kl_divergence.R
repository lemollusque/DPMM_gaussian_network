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

source("dao.R")
source("fns.R")

logdet <- function(S) {
  determinant(S, logarithm = TRUE)$modulus[1]
}

kl_cov <- function(S1, S2){
  p <- nrow(S1)
  0.5 * (
    sum(diag(solve(S2, S1))) -
      p +
      (logdet(S2) - logdet(S1))
  )
}

kl_sym <- function(S1, S2){
  0.5 * (kl_cov(S1, S2) + kl_cov(S2, S1))
}
init.seed <- 100
iter <- 3

param_grid <- expand.grid(
  lb_b = c(0, 0.5, 1), 
  ub_b = c(1, 1.5, 2), 
  lb_o = c(0.25, 0.5, 1), 
  ub_o = c(1, 1.5, 2), 
  lb_b2 = c(0, 0.5, 1), 
  ub_b2 = c(1, 1.5, 2), 
  lb_o2 = c(0.25, 0.5, 1), 
  ub_o2 = c(1, 1.5, 2), 
  n = 4
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
    .packages = c("matrixStats", "dplyr", "mclust")
  ) %dorng% {
    
    source("dao.R")
    source("fns.R")

    # set params
    j <- sim_grid$j[k]
    i <- sim_grid$i[k]
    
    lb_b <- param_grid$lb_b[j]
    ub_b <- param_grid$ub_b[j]
    lb_o <- param_grid$lb_o[j]
    ub_o <- param_grid$ub_o[j]
    
    lb_b2 <- param_grid$lb_b2[j]
    ub_b2 <- param_grid$ub_b2[j]
    lb_o2 <- param_grid$lb_o2[j]
    ub_o2 <- param_grid$ub_o2[j]
    
    n <- param_grid$n[j]
    
    g <- er_dag(n, d=0.2)
    g <- sf_out(g)
    truegraph <- randomize_graph(g)
    
    model1 <- cov(truegraph,lb_b = lb_b, ub_b = ub_b, lb_o = lb_o, ub_o = ub_o)
    model2 <- cov(truegraph,lb_b = lb_b2, ub_b = ub_b2, lb_o = lb_o2, ub_o = ub_o2)
    
    cov_diff <- norm(model1$S - model2$S, "F")
    kl <- kl_cov(model1$S, model2$S)
    kl_rev <- kl_cov(model2$S, model1$S)
    kl_sym_val <- kl_sym(model1$S, model2$S)
    
    iter_results <- data.frame(n = 4, 
                               lb_b = lb_b, ub_b = ub_b, lb_o = lb_o, ub_o = ub_o,
                               lb_b2 = lb_b2, ub_b2 = ub_b2, lb_o2 = lb_o2, ub_o2 = ub_o2,
                               cov_diff, kl, kl_rev, kl_sym_val)
    
    # show progress
    p()
    
    iter_results
  }
})
future::plan(sequential)



saveRDS(results, "Results/sims_kl.rds")
# results <- as.data.frame(readRDS("Results/Sims_Results_dp.rds"))


summary_params <- results %>%
  group_by(lb_b, ub_b, lb_o, ub_o,
           lb_b2, ub_b2, lb_o2, ub_o2) %>%
  summarise(
    mean_kl = mean(kl_sym_val),
    sd_kl = sd(kl_sym_val),
    mean_cov_diff = mean(cov_diff),
    .groups = "drop"
  )

best <- summary_params %>%
  arrange(desc(mean_kl)) %>%
  head(10)

worst <- summary_params %>%
  arrange(mean_kl) %>%
  head(10)

ggplot(summary_params,
       aes(abs(lb_b2 - lb_b), abs(ub_b2 - ub_b), fill = mean_kl)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_bw()

