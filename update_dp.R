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
iter <- 10
dual <- TRUE

# dirichlet params
alpha_prior <- c(2, 4)
initial_clusters <- 10
dp_iter <- 10
burnin <- 9
L <- 1

N = 100
n = 10
d = 1
bge.par = 0.01


# generate data
myDAG <- pcalg::randomDAG(n, prob = 0.2, lB = 1, uB = 2) 
trueDAG <- as(myDAG, "matrix")
truegraph <- 1*(trueDAG != 0)
data <- Fou_nldata(truegraph, N, lambda = d, noise.sd = 1, standardize = T)
if (is.null(colnames(data))) {
  colnames(data) <- paste0("v", seq_len(ncol(data)))
}
    
# start searchspace
alpha = 0.05
cor_mat <- cor(data)
startspace <- dual_pc(cor_mat, nrow(data), alpha = alpha, skeleton = T)


needed_score_sets <- function(startspace) {
  p <- ncol(startspace)
  out <- list()
  
  for (j in seq_len(p)) {
    base <- which(startspace[, j] == 1)
    extra <- setdiff(seq_len(p), c(j, base))

    out[[length(out) + 1]] <- list(
      child = j,
      parents = base
    )
    
    for (e in extra) {
      out[[length(out) + 1]] <- list(
        child = j,
        parents = sort(c(base, e))
      )
    }
  }
  out
}

needed_sets = needed_score_sets(startspace)

for (i in seq_along(needed_sets)) {
  child = needed_sets[[i]]$child
  parents = needed_sets[[i]]$parents
  dp_data = data[,c(child, parents)]
  if (length(parents) == 0){
    dp <-  DirichletProcessGaussian(dp_data)
  }
  else{
    dp <-  DirichletProcessMvnormal(dp_data, numInitialClusters = 10)
  }
  dp <- Fit(dp, dp_iter)
  
  Gamma_sample <- dp_membership_probs(dp, burnin, L)
  # add meta data in list
  Gamma_list <- add_membershipp(Gamma_list, 
                                Gamma_sample, 
                                child=child, 
                                parents=parents, 
                                active=FALSE)
}
dp_usrpar <- list(
  pctesttype = "bge",
  am = bge.par,
  dp_iter = dp_iter,
  dp_burnin = burnin,
  dp_n_sample = L
)
score <- scoreparameters("usr", data, usrpar = dp_usrpar)
searchspace <- iterativeMCMC(scorepar = score, startspace = startspace, hardlimit = 14, 
                             verbose = F, scoreout = TRUE, alphainit = 0.01)

