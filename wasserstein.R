library(BiDAG)
library(Bestie)
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
library(BayesFactor)
library(matrixStats)

library(data.table) # for last
library(DiagrammeR) # for making DAG plot
library(DiagrammeRsvg) ## for exporting svg for plotting to file
library(rsvg) ## for converting svg to png

source("toyDAGfunctionsSachs.R")

source("comparison_algs.R")
source("Fourier_fns.R")
source("dualPC.R")
source("dao.R")
source("fns.R")
insertSource("fns.R", package = "BiDAG")


set.seed(9)
N = 100
n = 10

nDAGs <- 500
bge.par = 0.01




# generate graph
myDAG <- pcalg::randomDAG(n, prob = 0.2, lB = 1, uB = 2) 
trueDAG <- as(myDAG, "matrix")
truegraph <- 1 * (trueDAG != 0)

sim <- simulate_bimodal(
  dag = t(truegraph),
  n = N,
  bimodal_sep = 5,
  return_model = TRUE
)

truegraph = t(sim$detectable_truegraph)

data <- sim$data
if (is.null(colnames(data))) {
  colnames(data) <- paste0("v", seq_len(ncol(data)))
}

# plot the true effects
B1 = sim$model1$B
B2 = sim$model2$B

# find effect
Eff1 <- t(solve(diag(ncol(B1)) - B1))
Eff2 <- t(solve(diag(ncol(B2)) - B2))

rescale_effects <- function(Eff, sds) {
  D <- diag(sds)
  D %*% Eff %*% solve(D)
}

Eff1 <- rescale_effects(Eff1, sim$raw_sd)
Eff2 <- rescale_effects(Eff2, sim$raw_sd)


trueEffectSamples <- sampleTrueEffects(
  Eff1 = Eff1,
  Eff2 = Eff2,
  n1 = sim$n1,
  n2 = sim$n2,
  labels = colnames(truegraph),
  nSamples = 1000
)


tol <- 1e-10
effectGraph <- 1 * (
  (abs(Eff1) > tol) |
    (abs(Eff2) > tol)
)
diag(effectGraph) <- 0

##############################################
# wasserstein wit dp
# dirichlet params
dp_usrpar <- list(
  pctesttype = "bge",
  am = bge.par,
  dp_prior = list(strength = 0.1, discount = 0),
  dp_mcmc = list(niter = 5000, nburn = 3000, model="LS"),
  dp_n_sample = 20,
  dp_fits = 1,
  dp_fitspace = "full"
)

DP.searchspace <- set.searchspace(
  data,
  "DP",
  usrpar = dp_usrpar
)

# prepare effects
trueSamples = makeTrueDAGSamples(
  truegraph = truegraph,
  nDAGs = nDAGs
)

alleffs = computeEffects_mem(
  sampledDAGs = trueSamples,
  scoreObject = DP.searchspace$score,
  DP = TRUE
)


res_dp <- wasserstein_effect_avg(
  effects4plot = alleffs,
  trueEffects = trueEffectSamples,
  effectGraph = effectGraph,
  sortlabs = 1:n
)

res_dp$avgW

#######################################################################
# bge
bge.searchspace <- set.searchspace(data, "bge", 0.01)
bge.alleffs = computeEffects_mem(sampledDAGs = trueSamples,
                   scoreObject = bge.searchspace$score,
                   DP=F) 

bge.res_dp <- wasserstein_effect_avg(
  effects4plot = bge.alleffs,
  trueEffects = trueEffectSamples,
  effectGraph = effectGraph,
  sortlabs = 1:n
)

bge.res_dp$avgW
