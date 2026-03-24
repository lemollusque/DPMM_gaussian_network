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
source("Fourier_fns.R")
source("dao.R")
source("fns.R")
insertSource("fns.R", package = "BiDAG")

init.seed <- 100
iter <- 100
dual <- TRUE

# dirichlet params
alpha_prior <- c(2, 4)
initial_clusters <- 10
dp_iter <- 200
dp_fits <- 2
burnin <- 190
L <- 10

N = 100
n = 4
d = 1
bge.par = 0.01


# generate data
g <- er_dag(n)
g <- sf_out(g)
g <- randomize_graph(g)
truegraph = t(g)
data <- Fou_nldata(truegraph, N, lambda = d, noise.sd = 1, standardize = T)
if (is.null(colnames(data))) {
  colnames(data) <- paste0("v", seq_len(ncol(data)))
}
    
# start searchspace
alpha = 0.05
cor_mat <- cor(data)
startspace <- dual_pc(cor_mat, nrow(data), alpha = alpha, skeleton = T)

# Find nodes with any incoming OR outgoing edges
idx <- which(rowSums(startspace) > 0 | colSums(startspace) > 0)
# Get node names
nodes <- rownames(startspace)[idx]

# prepare dirichlet gamma list
Gamma_list <- list()
for (f in seq_len(dp_fits)) {
  dp <- DirichletProcessMvnormal(data[, nodes], 
                                 alphaPriors = alpha_prior,
                                 numInitialClusters = initial_clusters)
  dp <- Fit(dp, dp_iter, progressBar = TRUE)
  
  Gamma_sample <- dp_membership_probs(dp, burnin, L)
  Gamma_list <- add_membershipp(
    Gamma_list,
    Gamma_sample,
    child = colnames(data)[1],
    parents = colnames(data)[2:ncol(data)],
    active = TRUE
  )
}

dp_usrpar <- list(
  pctesttype = "bge",
  am = bge.par,
  membershipp_list = Gamma_list,
  edgepf = 1
)