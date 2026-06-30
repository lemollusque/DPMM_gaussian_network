
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
source("dualPC.R")

library(data.table) # for last
library(DiagrammeR) # for making DAG plot
library(DiagrammeRsvg) ## for exporting svg for plotting to file
library(rsvg) ## for converting svg to png

inputData <- scale(data)


# dirichlet params
dp_usrpar <- list(
  pctesttype = "bge",
  am = 0.01,
  Imat = Imat,
  dp_prior = list(strength = 1, discount = 0),
  dp_mcmc = list(niter = 4000, nburn = 3000, model="LS"),
  dp_n_sample = 100,
  dp_fits = 5,
  dp_fitspace = "full",
  bgremove = TRUE
)


score <- scoreparameters("usr", inputData, usrpar = dp_usrpar)

myDAG <- pcalg::randomDAG(ncol(inputData)+ncol(Imat), prob = 0.2, lB = 1, uB = 2) 
trueDAG <- as(myDAG, "matrix")
truegraph <- 1*(trueDAG != 0)
truegraph[, 1:ncol(Imat)] <- 0

DAGscore(score, truegraph)

DPscoreDAG(score, truegraph)


