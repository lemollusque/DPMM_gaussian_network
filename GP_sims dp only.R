library(BiDAG)
library(matrixStats)
library(dirichletprocess)
library(dplyr)
library(ggplot2)
library(foreach)
library(doParallel)
library(parallelly)
library(mclust)

source("Fourier_fns.R")
source("comparison_algs.R")
source("dualPC.R")
source("dao.R")
source("fns.R")
# load functions script
insertSource("fns.R", package = "BiDAG")



init.seed <- 100
iter <- 50  # number of simulations
lambdas <- c(1)  # non-linearity; zero is linear
dual <- F    # use dualPC
n <- 4    # number of nodes
N <- 100     # number of samples
results <- data.frame()

# Parameters for ROC curves
bge.mus <- c(1)

computing_len = length(bge.mus) * iter * length(lambdas)
computing_done = 1

sink("logfile.txt")
for(lambda in lambdas) {
  print(paste("Lambda:", lambda))

  for(k in 1:length(bge.mus)) {
    print(paste("Parameter:", k))
    bge.par <- bge.mus[k]
    alpha_w  <- n + bge.par + 1      
    t <- bge.par * (alpha_w - n - 1) / (bge.par + 1)
    
    for (i in 1:iter) {
      print(paste("Iteration:", i))
      print(paste("computation: ", computing_done/computing_len))
      set.seed(init.seed+i)
      
      # Generate DAG & data
      g <- er_dag(n, ad =3)
      g <- sf_out(g)
      truegraph <- randomize_graph(g)
      
      
      model1 <- corr(g)
      model2 <- corr(g)
      
      df <- 3
      t_err <- function(n, var) {
        rt(n, df=df) * sqrt(var * (df-2)/df)
      }
      X1 <- simulate(model1$B, model1$O, N/2, err = t_err)
      X2 <- simulate(model2$B, model2$O, N/2, err = t_err)
      
      shift <- runif(ncol(X2), -2, 2)
      X2 <- sweep(X2, 2, shift, "+")
      
      data <- standardize(rbind(X1, X2))
      
      # Set initial search spaces
      DP.searchspace = set.searchspace(data, 
                                       dual, 
                                       "DP", 
                                       usrpar = list(pctesttype = "bge",
                                                     am = bge.par, 
                                                     aw = alpha_w, 
                                                     T0scale = t,
                                                     dp_iter = 100,
                                                     burnin = 50,
                                                     L = 10,
                                                     edgepf = 1
                                       )
      )
      for (p in 1:length(DP.searchspace$score$dp_scoreparam_list)){
        for ( s in 1:length(DP.searchspace$score$dp_scoreparam_list[[p]]$scores)){
          if (DP.searchspace$score$dp_scoreparam_list[[p]]$scores[[s]]$K != 1){
            cat(DP.searchspace$score$dp_scoreparam_list[[p]]$scores[[s]]$K)
            
          }
        }
      }
      
      
      bge.searchspace = set.searchspace(data, dual, "bge", bgepar = list(am = bge.par))
      
      # Bge score, partition
      bge.fit <- bge.partition.mcmc(bge.searchspace, order = F)
      results <- compare_results(bge.fit, c(bge.par, "BGe, partition", lambda), results, truegraph)
      
      # Bge score, order
      bge.fit <- bge.partition.mcmc(bge.searchspace, order = T)
      results <- compare_results(bge.fit, c(bge.par, "BGe, order", lambda), results, truegraph)
      
      # DP score, partition
      dp.fit <- DP.partition.mcmc(DP.searchspace, order = F)
      results <- compare_results(dp.fit, c(bge.par, "DP, partition", lambda), results, truegraph)
      
      # DP score, order
      dp.fit <- DP.partition.mcmc(DP.searchspace, order = T)
      results <- compare_results(dp.fit, c(bge.par, "DP, order", lambda), results, truegraph)
      
      # for printing
      computing_done = computing_done + 1
    }
  }
}
colnames(results) <- c("ESHD", "eTP", "eFP", "TPR", "FPR_P", 
                       "time", "parameter", "method", "lambda", "graph")
saveRDS(results, "Results/Sims_Results_dp.rds")
sink()


# plots resutts
keep_methods <- c("DP, partition", "DP, order", "BGe, partition", "BGe, order")

results_small <- results %>%
  filter(graph == "pattern", method %in% keep_methods) %>%
  mutate(
    lambda = factor(as.character(round(as.numeric(lambda), 2)),
                    levels = c("1", "0.5", "0")),
    method = factor(method, levels = keep_methods)
  )
results_small <- results_small %>%
  mutate(ESHD = as.numeric(ESHD))

lab_lambda <- function(x) parse(text = paste0("lambda==", x))
ggplot(results_small, aes(x = method, y = ESHD, color = method)) +
  geom_boxplot(aes(group = method), width = 0.6, outlier.shape = NA, linewidth = 0.6) +
  geom_jitter(aes(group = method), width = 0.15, alpha = 0.7, size = 1.2) +
  facet_wrap(~ lambda, nrow = 1, labeller = labeller(lambda = \(x) paste0("λ = ", x))) +
  coord_cartesian(ylim = c(0, 12)) +
  labs(x = NULL, y = "E=SHD") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major.x = element_blank()
  )
