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

init.seed <- 234
set.seed(init.seed)

bge.par = 0.01
# dirichlet params
dp_usrpar <- list(
  pctesttype = "bge",
  am = bge.par,
  dp_prior = list(strength = 0.0005, discount = 0),
  dp_mcmc = list(niter = 4000, nburn = 3000, model="LS"),
  dp_n_sample = 100,
  dp_fits = 1,
  dp_fitspace = "full"
)
output <- list(out_param = TRUE, out_type = "FULL")  


cor_mat <- cor(sachs.data)
fitspace <- dual_pc(cor_mat, nrow(sachs.data), alpha = 0.05, skeleton = T)
child = "PIP2"
parents <- names(which(fitspace[ , child] == 1))
children <- names(which(fitspace[child, ] == 1))
neighbour <- setdiff(unique(c(parents, children)), child)

dp_data = sachs.data[,c(child, neighbour)]
fit <- PYdensity(y = dp_data, mcmc = dp_usrpar$dp_mcmc, prior = dp_usrpar$dp_prior, output = output)

part <- partition(fit)
hard_clusters <- part$partitions[1, ]
table(hard_clusters)

data = sachs.data[hard_clusters == 1,]
dname = "cluster1"

nDAGs <- 50
nSeeds <- 50
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
  
  print(paste("Seed is", seednumber))
  
  
  load(file = paste0("./Sachs/saveout_subpopulation/dagdraw", nNodes, "seed", seednumber, "", ".RData"))
  
  bge.searchspace <- set.searchspace(data, "bge", bge.par)
  causalMats <- DAGintervention(sampledDAGs, bge.searchspace$score, sample=TRUE)
  
  save(causalMats,
       file=paste0("./Sachs/saveout_subpopulation/effects", nNodes, "seed", seednumber, dname, ".RData"))
  
  TRUE
}

nSeeds <- length(batch)
alldigraphs <- vector("list", nDAGs * nSeeds) # to store the graphs
alleffs <- vector("list", nDAGs * nSeeds) # to store the matrices of effects
for (nlevel in 1:nSeeds) {
  seednumber <- batch[nlevel]
  ## Retrieve sampled DAGs - DAG chain
  load(file = paste0("./Sachs/saveout_subpopulation/dagdraw", nNodes, "seed", seednumber, "", ".RData"))
  alldigraphs[1:nDAGs + (nlevel - 1) * nDAGs] <- sampledDAGs # remove the starting point
  ## Retrieve estimated effects
  load(file = paste0("./Sachs/saveout_subpopulation/effects", nNodes, "seed", seednumber, dname, ".RData"))
  alleffs[1:nDAGs + (nlevel - 1) * nDAGs] <- causalMats # remove the starting point
}
data4plot = list(alldigraphs=alldigraphs, alleffs=alleffs)

pdf(paste0("Sachs/SachsEffects_", dname, ".pdf"), width = 6, height = 6)
plotEffects(effects4plot = data4plot$alleffs, xmargs = c(0.1, 0.3), label_size = 1.5,
            sortlabs = 1:11, title_text = "")
dev.off()


