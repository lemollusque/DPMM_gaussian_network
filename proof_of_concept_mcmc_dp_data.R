## ---------------------------------------------------------------------------------------------------------------------------------------------------------
## load libraries
## ---------------------------------------------------------------------------------------------------------------------------------------------------------
require(BiDAG) ## for DAG sampling
require(Bestie) ## for intervention effects

require(dplyr)
require(tidyverse)
require(magrittr)
require(foreach)
require(doParallel)
require(parallelly)
require(DiagrammeR)
require(DiagrammeRsvg) ## for exporting svg for plotting to file
require(rsvg) ## for converting svg to png
require(dirichletprocess)

# load functions script
source("fns.R")


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
## Load and prepare the data
## ---------------------------------------------------------------------------------------------------------------------------------------------------------
dataname = "mixture"
set.seed(12)

gen_switching_sem <- function(N = 200, beta = 3, noise_sd = 0.3, standardize = TRUE, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # exogenous
  x2 <- rnorm(N)
  x3 <- rnorm(N)
  
  # x3 -> x1
  x1 <- 1.5 * x3 + rnorm(N, sd = noise_sd)
  
  # latent regime (+1 or -1)
  s <- sample(c(-1, 1), size = N, replace = TRUE)
  
  # x1 -> x4 (sign switches), and x2 -> x4
  x4 <- (s * beta) * x1 + 1.5 * x2 + rnorm(N, sd = noise_sd)
  
  dat <- data.frame(x1 = x1, x2 = x2, x3 = x3, x4 = x4)
  
  if (standardize) dat <- as.data.frame(scale(dat))
  attr(dat, "regime") <- s  # keep regimes if you want to inspect
  dat
}

n <- 4      # number of nodes
N <- 100    # number of samples
vars  <- c("x1","x2","x3","x4")
truegraph <- matrix(c(
  0,0,0,1,  # x1 -> x4
  0,0,0,1,  # x2 -> x4
  1,0,0,0,  # x3 -> x1
  0,0,0,0
), 4, byrow=TRUE, dimnames=list(vars, vars))

data <- gen_switching_sem(N = 200, beta = 3, noise_sd = 0.3, seed = 123, standardize = TRUE)

## ---------------------------------------------------------------------------------------------------------------------------------------------------------
## End data
## ---------------------------------------------------------------------------------------------------------------------------------------------------------

# Initiate params for DP and BGe
n <- ncol(data)
alpha_mu <- 1          
alpha_w  <- n + alpha_mu + 1      
t <- alpha_mu * (alpha_w - n - 1) / (alpha_mu + 1)

# DP
g0Priors <- list(
  mu0    = rep(0, n),
  Lambda = diag(n) / t,   # T = (1/t) I
  kappa0 = alpha_mu,
  nu     = alpha_w
)

scaled_data = scale(data) 
n_iter = 200
burnin = 180
L = 10 # sample to take

Gamma_list <- list()
vars  <- c("x1","x2","x3","x4")
for (child in vars){
  parents <- vars[vars != child]
  dp_data = scaled_data[,c(child, parents)]
  dp <-  DirichletProcessMvnormal(dp_data, g0Priors)
  dp <- Fit(dp, n_iter)
  
  Gamma_sample <- dp_membership_probs(dp, n_iter, burnin, L)
  Gamma_list <- add_membershipp(Gamma_list, Gamma_sample, child=child, parents=parents)
}


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
## DAG sampling
## ---------------------------------------------------------------------------------------------------------------------------------------------------------
nDAGs <- 100
batch <- 101:102
nSeeds <- length(batch)
labels4plot <- colnames(data)


# use parallelization for sampling
numCores <- min(length(batch), parallelly::availableCores())
cl <- makeCluster(numCores)

clusterEvalQ(cl, {
  library(BiDAG)
  #----------------------  overwrite functions ----------------------------------
  source("fns.R")  # must define usrscoreparameters + usrDAGcorescore replacements
  
  unlockBinding("usrscoreparameters", asNamespace("BiDAG"))
  assign("usrscoreparameters", usrscoreparameters, envir = asNamespace("BiDAG"))
  lockBinding("usrscoreparameters", asNamespace("BiDAG"))
  
  unlockBinding("usrDAGcorescore", asNamespace("BiDAG"))
  assign("usrDAGcorescore", usrDAGcorescore, envir = asNamespace("BiDAG"))
  lockBinding("usrDAGcorescore", asNamespace("BiDAG"))
})

# also make RNG reproducible across workers
parallel::clusterSetRNGStream(cl, 123)
registerDoParallel(cl)
# sampling loop
foreach(seednumber=batch) %dopar% {
  timing <- proc.time()
  print(paste("Seed is", seednumber))
  sampleDAGs(inData=scaled_data,
             scoretype="usr",
             usrpar = list(pctesttype = "bge",
                           membershipp_list = Gamma_list,
                           am = alpha_mu, 
                           aw = alpha_w, 
                           T0scale = t,
                           edgepf = 1
             ),
             nDigraphs=nDAGs,
             seed=seednumber,
             dname=dataname)
  print(proc.time() - timing)
}
stopCluster(cl)

## ---------------------------------------------------------------------------------------------------------------------------------------------------------
## Draw DAG
## ---------------------------------------------------------------------------------------------------------------------------------------------------------
alldigraphs <- vector("list", nDAGs * nSeeds) # to store the graphs
for (nlevel in 1:nSeeds) {
  seednumber <- batch[nlevel]
  ## Retrieve sampled DAGs - DAG chain
  load(file = paste0("./saveout/dagdraw", n, "seed", seednumber, dataname, ".RData"))
  alldigraphs[1:nDAGs + (nlevel - 1) * nDAGs] <- sampledDAGs[-1] # remove the starting point
}

# plot the DAG
graph2plot <- dagviz(alldigraphs,
                     style_mat = matrix(1, nrow=2*n + 1, ncol=2*n)[1:n, 1:n],
                     )
displayDAG(g2plot = graph2plot, figname=paste0("./plots/DAG_", dataname, ".png"))
