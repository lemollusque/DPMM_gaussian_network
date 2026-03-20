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

set.seed(1)

N = 1000
n = 4
g <- er_dag(n)
g <- sf_out(g)
truegraph <- randomize_graph(g)

model <- corr(truegraph)

X <- simulate_bimodal(model$B, model$O, n=N, bimodal_sep=6)
data <- standardize(X)

# ---- marginal densities ----
op <- par(mfrow=c(2, 2), mar=c(4, 4, 3, 1))

for (i in 1:ncol(data)) {
  hist(data[, i],
       probability=TRUE,
       breaks=30,
       main=paste("Distribution of X", i),
       xlab=paste("X", i),
       col="lightgray")
  lines(density(data[, i]), lwd=2)
}

par(op)

# ---- pairwise 4x4 plot ----
panel.density <- function(x, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5))
  d <- density(x)
  y <- d$y / max(d$y)
  lines(d$x, y, lwd=2)
}

panel.scatter <- function(x, y, ...) {
  points(x, y, pch=16, cex=0.5, ...)
}

pairs(data,
      diag.panel=panel.density,
      lower.panel=panel.scatter,
      upper.panel=panel.scatter,
      main="Pairwise distributions of the 4 variables")

print(t(truegraph))
