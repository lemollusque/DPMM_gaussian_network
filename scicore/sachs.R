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
  "data.table"
)

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

source("toyDAGfunctionsSachs.R")

source("comparison_algs.R")
source("Fourier_fns.R")
source("dualPC.R")
source("dao.R")
source("fns.R")
insertSource("fns.R", package = "BiDAG")

sachs.data <- read.csv("Sachs/2005_sachs_2_cd3cd28icam2_log_std.csv")
sachs.data <- as.matrix(sachs.data)

bge.par = 0.01
# dirichlet params
dp_usrpar <- list(
  pctesttype = "bge",
  am = bge.par,
  dp_prior = list(strength = 1, discount = 0),
  dp_mcmc = list(niter = 4000, nburn = 3000, model="LS"),
  dp_n_sample = 100,
  dp_fits = 1,
  dp_fitspace = "full"
)

nDAGs <- 50
nSeeds <- 50
batch <- 100 + 1:nSeeds
labels4plot <- colnames(sachs.data) 
nNodes <- length(labels4plot)

plan(multisession, workers = min(length(batch), availableCores()))
registerDoFuture()

foreach(
  seednumber = batch,
  .packages = c("BiDAG", "Bestie", "data.table", "mvtnorm")
) %dorng% {
  
  source("fns.R")
  insertSource("fns.R", package = "BiDAG")
  source("toyDAGfunctionsSachs.R")
  
  timing <- proc.time()
  print(paste("Seed is", seednumber))
  
  DP.searchspace <- set.searchspace(
    sachs.data,
    "DP",
    usrpar = dp_usrpar
  )
  
  sampleDAGs(
    inData = sachs.data,
    searchspace = DP.searchspace,
    weighted = TRUE,
    nDigraphs = nDAGs,
    seed = seednumber
  )
  
  computeEffects(
    n = nNodes,
    seed = seednumber,
    DP = TRUE
  )
  
  TRUE
}

