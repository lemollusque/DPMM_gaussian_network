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

# coefficients for chain X1 -> X2 -> X3
b21 <- 0.8
b32 <- 0.8

# error variances
o1 <- 1
o2 <- 1
o3 <- 1

# simulate chain
X1 <- bimodal_err(1000, o1, sep_sd=2)              # bimodal source
X2 <- b21 * X1 + rnorm(1000, 0, sqrt(o2))         # child of X1
X3 <- b32 * X2 + rnorm(1000, 0, sqrt(o3))         # child of X2

# plot
op <- par(mfrow=c(2, 3), mar=c(4, 4, 3, 1))

hist(X1, probability=TRUE, breaks=40,
     main="X1: source (bimodal)", xlab="X1", col="lightgray")
lines(density(X1), lwd=2)

hist(X2, probability=TRUE, breaks=40,
     main="X2: child of X1", xlab="X2", col="lightgray")
lines(density(X2), lwd=2)

hist(X3, probability=TRUE, breaks=40,
     main="X3: grandchild of X1", xlab="X3", col="lightgray")
lines(density(X3), lwd=2)

plot(X1, X2, pch=16, cex=0.4,
     main="Parent vs child: X1 -> X2",
     xlab="X1", ylab="X2")

plot(X2, X3, pch=16, cex=0.4,
     main="Parent vs child: X2 -> X3",
     xlab="X2", ylab="X3")

plot(X1, X3, pch=16, cex=0.4,
     main="Source vs grandchild: X1 -> X3",
     xlab="X1", ylab="X3")

par(op)


####################################graph
N = 1000
n = 4
g <- er_dag(n)
g <- sf_out(g)
truegraph <- randomize_graph(g)

model <- corr(truegraph)

out <- simulate_bimodal(model$B, model$O, n=N, bimodal_sep=2)
out$bimodal_node
X <- out$X
    
