# load libraries
packages <- c(
  "BiDAG",
  "dplyr",
  "ggplot2",
  "foreach",
  "doFuture",
  "future",
  "progressr",
  "doRNG",
  "mvtnorm",
  "BayesFactor",
  "matrixStats",
  "BNPmix",
  "data.table",
  "mclust",
  "pcalg",
  "graph"
)

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# preprocessing etc of the Sachs data
source("./isachssetup.R")


# Use BiDAG with intervention scoring
source("ifnsdp.R")
insertSource("ifnsdp.R", package = "BiDAG")

# load causal pipeline taken and adapted from https://github.com/annlia/causalpipe
source("itoyDAGfunctionsSachs.R")
source("intfns.R")

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
  dp_n_sample = 100,
  dp_fits = 2,
  dp_fitspace = "full",
  bgremove = TRUE
)

plan(multisession, workers = min(length(batch), availableCores()))
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
