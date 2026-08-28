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
source("comparison_algs.R")

library(data.table) # for last
library(DiagrammeR) # for making DAG plot
library(DiagrammeRsvg) ## for exporting svg for plotting to file
library(rsvg) ## for converting svg to png

inputData <- scale(data)
nDAGs <- 50
nSeeds <- 50
batch <- 100 + 1:nSeeds
labels4plot <- colnames(inputData) 
nNodes <- length(labels4plot)

dp_usrpar <- list(
  pctesttype = "bge",
  am = 0.01,
  Imat = Imat,
  dp_prior = list(strength = 0.0002, discount = 0),
  dp_mcmc = list(niter = 4000, nburn = 3000, model="LS"),
  dp_n_sample = 50,
  dp_fits = 2,
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
  source("itoyDAGfunctionsSachs.R")
  
  timing <- proc.time()
  print(paste("Seed is", seednumber))
  sampleDAGs(inData=inputData, scoretype = "usr",
             usrpar = dp_usrpar,
             nDigraphs=nDAGs, seed=seednumber,
             weighted=TRUE)
  computeEffects(n=nNodes, seed=seednumber, DP=TRUE)
  print(proc.time() - timing)
  TRUE
}

data4plot <- loadsamples(seeds=batch, nn=nNodes)

graph2plot <- dagviz(data4plot$alldigraphs, rm_nodes = 1:ncol(Imat), style_mat = matrix(1, 11, 11), title_text = "")
rsvg_png(charToRaw(export_svg(graph2plot)), "SachsDAGs.png", width = 4000)

pdf("SachsEffects.pdf", width = 6, height = 6)
plotEffects(effects4plot = data4plot$alleffs, xmargs = c(0.1, 0.3), label_size = 1.5,
            sortlabs = 1:11, title_text = "")
dev.off()

# edge proba
edge_probabilities <- function(dags) {
  prob_mat <- Reduce(
    "+",
    lapply(dags, function(x) as.matrix(x) * 1)
  ) / length(dags)
  
  diag(prob_mat) <- NA_real_
  
  as.data.frame(
    as.table(prob_mat),
    responseName = "probability"
  ) %>%
    filter(!is.na(probability)) %>%
    transmute(
      from = Var1,
      to = Var2,
      probability
    )
}

bgn <- ncol(Imat)
sampledDAGs <- lapply(data4plot$alldigraphs, function(A) {
  A[-(1:bgn), -(1:bgn), drop = FALSE]
})


edge_probs <- edge_probabilities(sampledDAGs)
edge_summary <- edge_probs %>%
  summarise(
    displayed_edges = sum(probability > 0.1),
    
    low_probability = sum(
      probability > 0.1 & probability < 0.5
    ),
    
    moderate_probability = sum(
      probability >= 0.5 & probability < 0.8
    ),
    
    high_probability = sum(
      probability >= 0.8
    ),
    
    very_high_probability = sum(
      probability >= 0.9
    ),
    
    proportion_high = high_probability / displayed_edges,
    proportion_very_high = very_high_probability / displayed_edges
  )

edge_summary
