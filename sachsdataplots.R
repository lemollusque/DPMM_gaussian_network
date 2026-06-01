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

trueDAG <- read.csv("Sachs/sachs_dag_with_missing.csv")
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
  upper = list(continuous = wrap("cor", size = 3, digits = 2)),
  lower = list(continuous = my_smooth),
  diag = list(continuous = wrap("barDiag", bins = 30, alpha = 0.5))
) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 7),
    strip.text = element_text(size = 8)
  )


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
truecpdag = pcalg::dag2cpdag(trueDAG)
truecpdag <- 1 * (truecpdag)
labels4plot <- colnames(sachs.data) 
nNodes <- length(labels4plot)
bge.searchspace <- set.searchspace(sachs.data, TRUE, "bge", 0.01)
nDAGs <- 50
makeTrueDAGSamples( truegraph = trueDAG, 
                    nDAGs = nDAGs, seed = 101, 
                    scoreObject = bge.searchspace$score, 
                    outdir = "./Sachs/saveout" ) 

computeEffects(
  n = nNodes,
  seed = 101
)
data4plot <- loadsamples(seeds=c(101), nn=nNodes)


graph2plot <- dagviz(data4plot$alldigraphs, style_mat = matrix(1, 11, 11), title_text = "")
rsvg_png(charToRaw(export_svg(graph2plot)), "Sachs/SachstrueDAG.png", width = 4000)
