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
alpha_prior = c(2, 4)
dp_iter = 100
dp_fits = 1
burnin = 90
L = 10

param_grid <- expand.grid(
  N = c(200),
  n = 4,
  d = c(1,2,5,10),
  bge.par = 1
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
    .packages = c("BiDAG", "matrixStats", "dirichletprocess", "dplyr", "mclust")
  ) %dorng% {
    
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
    
    g <- er_dag(n, d=0.2)
    g <- sf_out(g)
    truegraph <- randomize_graph(g)
    
    model1 <- cov(truegraph)
    model2 <- cov(truegraph)
    
    
    X1 <- rmvt(N/2, sigma =  model1$S, df = 3)
    X2 <- rmvt(N/2, sigma =  model2$S, df = 3)
    
    v <- rnorm(ncol(X2))
    v <- v / sqrt(sum(v^2))   
    shift <- d * v
    X2 <- sweep(X2, 2, shift, "+")
    
    data <- standardize(rbind(X1, X2))
    if(is.null(colnames(data))) {
      colnames(data) <- sapply(c(1:ncol(data)), function(x) paste("v", x, sep = ""))
    }
    
    # prepare dirichlet gamm list
    Gamma_list <- list()
    for (f in seq_len(dp_fits)){
      dp <-  DirichletProcessMvnormal(data, alphaPriors = alpha_prior)
      dp <- Fit(dp, dp_iter, progressBar = FALSE)
      
      Gamma_sample <- dp_membership_probs(dp, burnin, L)
      Gamma_list <- add_membershipp(Gamma_list, 
                                    Gamma_sample, 
                                    child=colnames(data)[1], 
                                    parents=colnames(data)[2:ncol(data)], 
                                    active=TRUE)
    }
    dp_usrpar = list(
      pctesttype = "bge",
      am = bge.par,
      membershipp_list = Gamma_list,
      edgepf = 1
    )
    
    # search spaces
    DP.searchspace <- set.searchspace(data, dual, "DP", usrpar = dp_usrpar)
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
    
    dp.fit <- DP.partition.mcmc(DP.searchspace, order = FALSE)
    iter_results <- compare_results(
      dp.fit, c(bge.par, "DP, partition"), iter_results, truegraph
    )
    
    dp.fit <- DP.partition.mcmc(DP.searchspace, order = TRUE)
    iter_results <- compare_results(
      dp.fit, c(bge.par, "DP, order"), iter_results, truegraph
    )
    
    iter_results$N <- N
    iter_results$n <- n
    iter_results$d <- d
    
    # show progress
    p()
    
    iter_results
  }
})
future::plan(sequential)

colnames(results) <- c(
  "ESHD", "eTP", "eFP", "TPR", "FPR_P",
  "time", "parameter", "method", "graph", "N", "n", "d"
)


saveRDS(results, "Results/Sims_Results_dp.rds")


# plots resutts
keep_methods <- c("DP, partition", "DP, order", "BGe, partition", "BGe, order")

results_small <- results %>%
  filter(graph == "pattern", method %in% keep_methods) 

results_small <- results_small %>%
  mutate(ESHD = as.numeric(ESHD))

# extract benchmark
data_benchmark <- as.data.frame(readRDS("Results/Sims_benchmark.rds"))
data_benchmark <- data_benchmark %>%
  filter(graph == "pattern", 
         method %in% keep_methods,
         n==param_grid$n[1],
         N==param_grid$N[1]) 
data_benchmark <- data_benchmark %>%
  mutate(ESHD = as.numeric(ESHD))

# add single multivariate gaussian benchmark
benchmark_medians <- data_benchmark %>%
  group_by(method) %>%
  summarise(median_ESHD = median(ESHD, na.rm = TRUE), .groups = "drop")

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
  geom_hline(
    data = benchmark_medians,
    aes(yintercept = median_ESHD, color = method),
    linetype = "dashed",
    linewidth = 1,
    show.legend = FALSE
  ) +
  
  coord_cartesian(ylim = c(0, 6)) +
  labs(x = NULL, y = "E=SHD") +
  facet_grid( ~ d, scales = "free_y") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major.x = element_blank()
  )

