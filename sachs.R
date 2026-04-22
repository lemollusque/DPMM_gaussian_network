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
library(readxl)

source("comparison_algs.R")
source("dualPC.R")
source("Fourier_fns.R")
source("dao.R")
source("fns.R")
insertSource("fns.R", package = "BiDAG")

sachs.data <- read_excel("Sachs/1. cd3cd28.xls")
trueDAGbn <- readRDS("Sachs/sachs_graph.rds")
set.seed(100)
sachs.data <- scale(log2(sachs.data + 0.5))
sachs.data <- standardize(sachs.data)

results <- data.frame()

bge.par = 0.01
dual <- TRUE
# dirichlet params
dp_iter <- 1000
burnin <- 500
L <- 100
dp_fits <- 4

# dp settings
dp_usrpar <- list(
  pctesttype = "bge",
  am = bge.par,
  dp_iter = dp_iter,
  dp_burnin = burnin,
  dp_n_sample = L,
  dp_fits = dp_fits
)

# search spaces
DP.searchspace <- set.searchspace(sachs.data, dual, "DP", usrpar = dp_usrpar)
saveRDS(DP.searchspace, "Sachs/dpsearchspace.rds")
DP.searchspace <- readRDS("Sachs/dpsearchspace.rds")

bge.searchspace = set.searchspace(sachs.data, dual, "bge", bge.par)

# DP score, partition
dp.fit <- DP.partition.mcmc(DP.searchspace, order = FALSE, iterations = 1200)
dp.edgep <- post.edges(dp.fit)
results <- compare_results(dp.fit, c(dp.edgep, "DP, partition"), results, trueDAGbn)

# DP score, order
dp.fit <- DP.partition.mcmc(DP.searchspace, order = TRUE, iterations = 1200)
dp.edgep <- post.edges(dp.fit)
results <- compare_results(dp.fit, c(dp.edgep, "DP, partition"), results, trueDAGbn)

# BGe score, partition
bge.fit <-  bge.partition.mcmc(bge.searchspace, order = F, iterations = 1200)
bge.edgep <- post.edges(bge.fit)
results <- compare_results(bge.fit, c(bge.edgep, "BGe, partition"), results, trueDAGbn)

# BGe score, order
bge.fit <-  bge.partition.mcmc(bge.searchspace, order = T, iterations = 1200)
bge.edgep <- post.edges(bge.fit)
results <- compare_results(bge.fit, c(bge.edgep, "BGe, order"), results, trueDAGbn)


colnames(results) <- c("ESHD", "eTP", "eFP", "TPR", "FPR_P", "time",
                       "ErktoAkt", "ErktoPKA", "Scorefn", "graph")
saveRDS(results, "Sachs/Sachs_results.rds")
