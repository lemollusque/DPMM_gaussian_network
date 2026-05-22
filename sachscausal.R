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

sachs.data <- read.csv("Sachs/2005_sachs_2_cd3cd28icam2_log_std.csv")
sachs.data <- as.matrix(sachs.data)

nDAGs <- 100
nSeeds <- 20
batch <- 100 + 1:nSeeds
labels4plot <- colnames(sachs.data) 
nNodes <- length(labels4plot)


plan(multisession, workers = min(length(batch), availableCores() - 1))
registerDoFuture()

foreach(
  seednumber = batch,
  .packages = c("BiDAG", "Bestie", "data.table", "mvtnorm")
) %dorng% {
  
  source("fns.R")
  insertSource("fns.R", package = "BiDAG")
  source("toyDAGfunctionsSachs.R")
  
  timing <- proc.time()
  
  fname <- sprintf(
    "Sachs/parallel_searchspaces/DP_searchspace_rep_%03d.rds",
    seednumber - 100
  )
  
  cat("Seed is", seednumber, "\n")
  
  DP.searchspace <- readRDS(fname)
  
  sampleDAGs(
    inData = sachs.data,
    scoreObject = DP.searchspace$score,
    nDigraphs = nDAGs,
    seed = seednumber
  )
  
  computeEffects(
    n = nNodes,
    seed = seednumber,
    DP = TRUE
  )
  
  cat("Time:", proc.time() - timing, "\n")
  
  TRUE
}

data4plot <- loadsamples(seeds=batch, nn=nNodes)

graph2plot <- dagviz(data4plot$alldigraphs, style_mat = matrix(1, 11, 11), title_text = "")
rsvg_png(charToRaw(export_svg(graph2plot)), "Sachs/SachsDAGs.png", width = 4000)

pdf("Sachs/SachsEffects.pdf", width = 6, height = 6)
plotEffects(effects4plot = data4plot$alleffs, xmargs = c(0.1, 0.3), label_size = 1.5,
            sortlabs = 1:11, title_text = "")
dev.off()


