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
library(transport)

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
# plot
pdf("effects/TrueEffects.pdf", width = 6, height = 6)
plotEffects(
  effects4plot = trueEffectSamples, xmargs = c(0.1, 0.3), label_size = 1.5,
  sortlabs = 1:n,
  title_text = ""
)
dev.off()


truegraph
Eff1
Eff2 
sim$n1
sim$n2


##############################################
# with perfect dags
##############################################
nDAGs <- 500

bge.par = 0.01
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

# plot effects
pdf("effects/EstimatedEffects.pdf", width = 6, height = 6)
plotEffects(effects4plot = alleffs, xmargs = c(0.1, 0.3), label_size = 1.5,
            sortlabs = 1:n, 
            title_text = "")
dev.off()

pdf("effects/Wasserstein.pdf", width = 6, height = 6)
# compare with true effects
Wdist <- plotCompareEffects(
  effects4plot = alleffs,
  trueEffects = trueEffectSamples,
  sortlabs = 1:n,
  title_text = ""
)
dev.off()

#######################################################################
# bge
bge.searchspace <- set.searchspace(data, "bge", 0.01)
bge.alleffs = computeEffects_mem(sampledDAGs = trueSamples,
                   scoreObject = bge.searchspace$score,
                   DP=F) 

pdf("effects/BGeEstimatedEffects.pdf", width = 6, height = 6)
plotEffects(effects4plot = bge.alleffs, xmargs = c(0.1, 0.3), label_size = 1.5,
            sortlabs = 1:n, 
            title_text = "")
dev.off()

pdf("effects/BGeWasserstein.pdf", width = 6, height = 6)
# compare with true effects
Wdist <- plotCompareEffects(
  effects4plot = bge.alleffs,
  trueEffects = trueEffectSamples,
  sortlabs = 1:n,
  title_text = ""
)
dev.off()
