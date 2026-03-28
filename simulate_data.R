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
source("Fourier_fns.R")
source("dualPC.R")
source("dao.R")
source("fns.R")
insertSource("fns.R", package = "BiDAG")

set.seed(1)

N = 100
n = 10


plot_simulated_data <- function(data, title){
  # ---- pairwise 4x4 plot ----
  panel.density <- function(x, ...) {
    usr <- par("usr")
    on.exit(par(usr = usr))
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
        main=title)
}

# generate graph
myDAG <- pcalg::randomDAG(n, prob = 0.2, lB = 1, uB = 2) 
trueDAG <- as(myDAG, "matrix")
truegraph <- 1 * (trueDAG != 0)

for(row in 1:ncol(truegraph)){
  cat(truegraph[row,], "\n")
}
# run bimodal data
data <- simulate_bimodal(t(truegraph), n=N, bimodal_sep=2)
plot_simulated_data(data, "Bimodal data")

# run bimodal on single node data
data <- simulate_bimodal_one_node(t(truegraph), n=N, bimodal_sep=2)
plot_simulated_data(data, "One node Bimodal data")

# fourier data
data <- Fou_nldata(truegraph, N, lambda = 1, noise.sd = 1, standardize = T, concentration = 1)
plot_simulated_data(data, "Fourier data")
