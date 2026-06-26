# preprocessing etc of the Sachs data
source("./isachssetup.R")

# load libraries
library(BiDAG)
library(Bestie)

library(matrixStats)
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
library(BNPmix)

# Use BiDAG with intervention scoring
source("ifnsdp.R")
insertSource("ifnsdp.R", package = "BiDAG")

# load causal pipeline taken and adapted from https://github.com/annlia/causalpipe
source("itoyDAGfunctionsSachs.R")
source("intfns.R")

library(data.table) # for last
library(DiagrammeR) # for making DAG plot
library(DiagrammeRsvg) ## for exporting svg for plotting to file
library(rsvg) ## for converting svg to png

inputData <- scale(data)
nDAGs <- 50
nSeeds <- 10
batch <- 100 + 1:nSeeds
labels4plot <- colnames(inputData) 
nNodes <- length(labels4plot)

dp_usrpar <- list(
  pctesttype = "bge",
  am = 0.01,
  Imat = Imat,
  dp_prior = list(strength = 0.0002, discount = 0),
  dp_mcmc = list(niter = 500, nburn = 300, model="LS"),
  dp_n_sample = 5,
  dp_fits = 1,
  dp_fitspace = "full",
  bgremove = TRUE
)

plan(multisession, workers = min(length(batch), availableCores()-1))
registerDoFuture()

foreach(
  seednumber = batch,
  .packages = c("BiDAG", "Bestie", "data.table", "mvtnorm", "BNPmix")
) %dorng% {
  
  source("ifnsdp.R")
  insertSource("ifnsdp.R", package = "BiDAG")
  source("intfns.R")
  
  timing <- proc.time()
  print(paste("Seed is", seednumber))
  sampleDAGs(inData=inputData, scoretype = "usr",
             usrpar = dp_usrpar,
             nDigraphs=nDAGs, seed=seednumber,
             weighted=FALSE)
  computeEffects(n=nNodes, seed=seednumber, DP=TRUE)
  print(proc.time() - timing)
  TRUE
}

data4plot <- loadsamples(seeds=batch, nn=nNodes)

graph2plot <- dagviz(data4plot$alldigraphs, rm_nodes = 1:6, style_mat = matrix(1, 11, 11), title_text = "")
rsvg_png(charToRaw(export_svg(graph2plot)), "SachsDAGs.png", width = 4000)

pdf("SachsEffects.pdf", width = 6, height = 6)
plotEffects(effects4plot = data4plot$alleffs, xmargs = c(0.1, 0.3), label_size = 1.5,
            sortlabs = 1:11, title_text = "")
dev.off()


