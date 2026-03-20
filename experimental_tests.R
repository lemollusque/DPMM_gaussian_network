library(BiDAG)
library(matrixStats)
library(dirichletprocess)
library(dplyr)
library(ggplot2)
library(foreach)
library(doFuture)
library(future)
library(parallelly)
library(aricode)
library(mclust)
library(openxlsx)
library(progressr)
library(mvtnorm)
library(doFuture)
library(doRNG)

source("dao.R")
source("dualPC.R")
source("fns.R")




N = 1000
shift_size=1
g <- er_dag(n)
g <- sf_out(g)
truegraph <- randomize_graph(g)

model <- corr(truegraph)
data_simulation <- simulate_bimodal_for_tuning(model$B, model$O, n=N, bimodal_sep=shift_size)
data <- standardize(data_simulation$X)
colnames(data) <- paste0("v", seq_len(ncol(data)))

truth <- ifelse(data_simulation$cluster == -1, "X1", "X2")




start <- Sys.time()
dp <- DirichletProcessMvnormal(
  data,
  numInitialClusters = 10
)

dp <- Fit(dp, dp_iter)
time <- Sys.time() - start
print(time)

print(dp$pointsPerCluster)
View(dp$weightsChain)




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


##################################

library(dirichletprocess)

set.seed(1)

# -----------------------------
# 1) Simulate simple bimodal data
# -----------------------------
N <- 100

z <- sample(1:2, N, replace = TRUE)
x1 <- rnorm(N, mean = ifelse(z == 1, -2, 2), sd = 0.8)  # bimodal variable
x2 <- rnorm(N)
x3 <- rnorm(N)
x4 <- rnorm(N)

data <- scale(cbind(x1, x2, x3, x4))
colnames(data) <- paste0("V", 1:4)

# -----------------------------
# 2) Same prior except kappa0
# -----------------------------
d <- ncol(data)
alpha_prior <- c(2, 2)

g0_large_kappa <- list(
  mu0 = colMeans(data),
  kappa0 = 4,          # stronger shrinkage
  nu = d + 2,
  Lambda = diag(d)
)

g0_small_kappa <- list(
  mu0 = colMeans(data),
  kappa0 = 0.05,       # much weaker shrinkage
  nu = d + 2,
  Lambda = diag(d)
)

# -----------------------------
# 3) Fit DP mixtures
# -----------------------------
dp_large <- DirichletProcessMvnormal(
  data,
  g0Priors = g0_large_kappa,
  alphaPriors = alpha_prior
)
dp_large <- Fit(dp_large, 200)

dp_small <- DirichletProcessMvnormal(
  data,
  g0Priors = g0_small_kappa,
  alphaPriors = alpha_prior
)
dp_small <- Fit(dp_small, 200)

# -----------------------------
# 4) Extract cluster labels
# -----------------------------
cl_large <- dp_large$clusterLabels
cl_small <- dp_small$clusterLabels

# If that fails in your version, run:
# names(dp_large)
# str(dp_large, max.level = 1)

# -----------------------------
# 5) Compare cluster assignments
# -----------------------------
op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

plot(data[,1], jitter(cl_large),
     col = cl_large, pch = 16, cex = 0.7,
     xlab = "Standardized V1",
     ylab = "Cluster label",
     main = "Large kappa0 = 4")

plot(data[,1], jitter(cl_small),
     col = cl_small, pch = 16, cex = 0.7,
     xlab = "Standardized V1",
     ylab = "Cluster label",
     main = "Small kappa0 = 0.05")

par(op)

# -----------------------------
# 6) Optional: show densities by true component
# -----------------------------
op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

plot(density(data[z == 1, 1]), lwd = 2,
     main = "True bimodality in V1", xlab = "V1")
lines(density(data[z == 2, 1]), lwd = 2, lty = 2)
legend("topright", legend = c("true comp 1", "true comp 2"),
       lty = c(1, 2), bty = "n")

plot(density(data[,1]), lwd = 2,
     main = "Marginal density of V1", xlab = "V1")

par(op)
