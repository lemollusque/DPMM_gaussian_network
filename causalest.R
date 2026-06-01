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

set.seed(7)
N = 500
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


sampleTrueEffects <- function(Eff1, Eff2, n1, n2, labels = NULL, nSamples = 1000) {
  n <- nrow(Eff1)
  ntot <- n1 + n2
  out <- vector("list", nSamples)
  for (s in seq_len(nSamples)) {
    M <- if (runif(1) < n1 / ntot) Eff1 else Eff2
    
    if (!is.null(labels)) {
      colnames(M) <- rownames(M) <- labels
    }
    out[[s]] <- M
  }
  out
}
trueSamples <- sampleTrueEffects(
  Eff1 = Eff1,
  Eff2 = Eff2,
  n1 = sim$n1,
  n2 = sim$n2,
  labels = colnames(sim$data),
  nSamples = 1000
)
# plot
pdf("TrueEffects.pdf", width = 6, height = 6)
plotEffects(
  effects4plot = trueSamples, xmargs = c(0.1, 0.3), label_size = 1.5,
  sortlabs = 1:n,
  title_text = "True Mixture Distribution of Causal Effects\n"
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
makeTrueDAGSamples <- function(truegraph, nDAGs, seed = 101, dname = "",
                               scoreObject, outdir = "./saveout") {
  
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  
  n <- ncol(truegraph)
  
  sampledDAGs <- replicate(
    nDAGs + 1,
    truegraph,
    simplify = FALSE
  )
  
  sampledDAGs <- lapply(sampledDAGs, function(A) {
    A <- as.matrix(A)
    colnames(A) <- rownames(A) <- colnames(scoreObject$data)
    A
  })
  
  DAGscores <- rep(NA_real_, length(sampledDAGs))
  
  save(
    sampledDAGs,
    DAGscores,
    scoreObject,
    file = file.path(outdir, paste0("dagdraw", n, "seed", seed, dname, ".RData"))
  )
}

nDAGs <- 50
nSeeds <- 50
seeds_per_dp <- 10
batch <- 100 + 1:nSeeds

dp_groups <- split(
  batch,
  ceiling(seq_along(batch) / seeds_per_dp)
)

labels4plot <- colnames(data) 
nNodes <- length(labels4plot)

# dp settings
dp_usrpar <- list(
  pctesttype = "bge",
  am = 0.01,
  dp_iter = 200,
  dp_burnin = 100,
  dp_n_sample = 10,
  dp_fits = 1,
  alphaPriors = c(8,4),
  g0Priors = function(n) {
    list(mu0 = rep(0, n), 
         kappa0 = 0.1, 
         nu = n+5, 
         Lambda = diag(n)*0.5)
  },
  numInitialClusters = min(20, ceiling(sqrt(N))),
  progressBar = T
)

plan(multisession, workers = min(length(dp_groups), parallelly::availableCores() - 1))
registerDoFuture()

handlers(global = TRUE)

foreach(
  g = seq_along(dp_groups),
  .packages = c("BiDAG", "matrixStats", "dirichletprocess", "dplyr", "mclust", "mvtnorm")
  ) %dorng% {
  
  source("toyDAGfunctionsSachs.R")
  source("comparison_algs.R")
  source("dualPC.R")
  source("Fourier_fns.R")
  source("dao.R")
  source("fns.R")
  insertSource("fns.R", package = "BiDAG")
  
  cat("Fitting DP group", g, "\n")
  
  DP.searchspace <- set.searchspace(
    data,
    "dual",
    "DP",
    usrpar = dp_usrpar
  )
  
  for (seednumber in dp_groups[[g]]) {
    
    cat("Seed is", seednumber, "\n")
    
    makeTrueDAGSamples(
      truegraph = truegraph,
      nDAGs = nDAGs,
      seed = seednumber,
      scoreObject = DP.searchspace$score,
      outdir = "./Sachs/saveout"
    )
    
    computeEffects(
      n = nNodes,
      seed = seednumber,
      DP = TRUE
    )
  }
  
  TRUE
}
data4plot <- loadsamples(seeds=batch, nn=nNodes)

pdf("EstimatedEffects.pdf", width = 6, height = 6)
plotEffects(effects4plot = data4plot$alleffs, xmargs = c(0.1, 0.3), label_size = 1.5,
            sortlabs = 1:n, 
            title_text = "Estimated Mixture Distribution of Causal Effects\n")
dev.off()

pdf("Wasserstein.pdf", width = 6, height = 6)
# compare with true effects
Wdist <- plotCompareEffects(
  effects4plot = data4plot$alleffs,
  trueEffects = trueSamples,
  sortlabs = 1:n,
  title_text = "Estimated Effects with Entrywise Wasserstein Distance\n"
)
dev.off()

#######################################################################
# bge

for (seednumber in batch) { 
  print(paste("Seed is", seednumber)) 
  # search spaces 
  bge.searchspace <- set.searchspace(data, TRUE, "bge", 0.01)
  makeTrueDAGSamples( truegraph = truegraph, 
                      nDAGs = nDAGs, seed = seednumber, 
                      scoreObject = bge.searchspace$score, 
                      outdir = "./Sachs/saveout" ) 
  computeEffects(n = nNodes, seed = seednumber, DP=F) 
}
data4plot <- loadsamples(seeds=batch, nn=nNodes)

pdf("BGeEstimatedEffects.pdf", width = 6, height = 6)
plotEffects(effects4plot = data4plot$alleffs, xmargs = c(0.1, 0.3), label_size = 1.5,
            sortlabs = 1:n, 
            title_text = "BGe Estimated Distribution of Causal Effects\n")
dev.off()

pdf("BGeWasserstein.pdf", width = 6, height = 6)
# compare with true effects
Wdist <- plotCompareEffects(
  effects4plot = data4plot$alleffs,
  trueEffects = trueSamples,
  sortlabs = 1:n,
  title_text = "BGe Estimated Effects with Entrywise Wasserstein Distance\n"
)
dev.off()
