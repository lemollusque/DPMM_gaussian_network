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
source("dao.R")
source("fns.R")
insertSource("fns.R", package = "BiDAG")

sachs.data <- read.csv("Sachs/2005_sachs_2_cd3cd28icam2_log_std.csv")

sachs.data <- standardize(sachs.data)
sachs.data <- as.matrix(sachs.data)
trueDAGbn <- read.csv("Sachs/sachs.csv")
trueDAGbn <- as.matrix(trueDAGbn)
set.seed(100)

results <- data.frame()

bge.par = 0.01
dual <- TRUE
# dirichlet params
dp_iter <- 1000
burnin <- 500
L <- 100
dp_fits <- 2

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
saveRDS(DP.searchspace, "Sachs/dpsearchspace_sachs_2.rds")
