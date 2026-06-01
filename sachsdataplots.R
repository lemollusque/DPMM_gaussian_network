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
library(GGally)
library(igraph)

source("comparison_algs.R")
source("dualPC.R")
source("dao.R")
source("fns.R")
insertSource("fns.R", package = "BiDAG")

sachs.data <- read.csv("Sachs/2005_sachs_2_cd3cd28icam2_log_std.csv")
sachs.data <- as.matrix(sachs.data)

trueDAG <- read.csv("Sachs/sachs.csv")
trueDAG <- as.matrix(trueDAG)

# plot data
my_smooth <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_point(color = "black", alpha = 0.5, size = 0.6) +
    geom_smooth(
      method = "lm",
      se = FALSE,
      color = "red",
      linewidth = 0.8
    )
}

ggpairs(
  sachs.data,
  upper = list(continuous = wrap("cor", size = 4)),
  lower = list(continuous = my_smooth),
  diag = list(continuous = wrap("barDiag", bins = 30, alpha = 0.5))
) +
  theme_minimal()

trueDAG <- read.csv("Sachs/sachs.csv")
trueDAG <- as.matrix(trueDAG)

# plot DAG
g <- graph_from_adjacency_matrix(
  trueDAG,
  mode = "undirected",
  diag = FALSE
)

plot(
  g,
  layout = layout_with_fr(g),   # Fruchterman-Reingold
  vertex.size = 25,
  vertex.color = "white",
  vertex.frame.color = "black",
  vertex.label.cex = 0.8,
  edge.width = 1.5
)