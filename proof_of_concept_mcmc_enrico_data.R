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
source("Fourier_fns.R")
source("dualPC.R")


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
## Load and prepare the data
## ---------------------------------------------------------------------------------------------------------------------------------------------------------
dataname = "fourier"
set.seed(101)
lambda <- 10
n <- 10  # number of nodes
N <- 100  # number of samples

myDAG <- pcalg::randomDAG(n, prob = 0.5, lB = 1, uB = 2) 
trueDAG <- as(myDAG, "matrix")
truegraph <- 1*(trueDAG != 0)
data <- Fou_nldata(truegraph, N, lambda = lambda, noise.sd = 1, standardize = T) 

vars <- paste0("x", 1:n)
colnames(data) = vars

# Initiate params for DP and BGe
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
burnin = 30
L = 10 # sample to take

cormat <- cor(scaled_data)
startspace <- dual_pc(cormat, nrow(scaled_data), alpha = 0.05, skeleton = T)

# start from fully connected:
startspace <- matrix(1, n, n)
diag(startspace) <- 0
dimnames(startspace) <- list(vars, vars)


Gamma_list <- list()
for (child in vars){
  parents <- names(which(startspace[ , child] == 1))
  dp_data = scaled_data[,c(child, parents)]
  if (length(parents) == 0){
    dp <-  DirichletProcessGaussian(dp_data,
                                    alphaPriors = c(1, 20))
  }
  else{
    n_col = ncol(dp_data)
    g0Priors <- list(
      mu0    = rep(0, n_col),
      Lambda = diag(n_col) / t,   # T = (1/t) I
      kappa0 = alpha_mu,
      nu     = alpha_w
    )
    
    dp <-  DirichletProcessMvnormal(dp_data, g0Priors)
  }
  dp <- Fit(dp, n_iter)
  
  Gamma_sample <- dp_membership_probs(dp, burnin, L)
  Gamma_list <- add_membershipp(Gamma_list, 
                                Gamma_sample, 
                                child=child, 
                                parents=parents, 
                                active=TRUE)
}
## ---------------------------------------------------------------------------------------------------------------------------------------------------------
## DAG sampling
## ---------------------------------------------------------------------------------------------------------------------------------------------------------
nDAGs <- 100
batch <- 101:120
nSeeds <- length(batch)
labels4plot <- colnames(data)


# use parallelization for sampling
numCores <- min(length(batch), parallelly::availableCores())
cl <- makeCluster(numCores)

clusterEvalQ(cl, {
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
  cat(paste("Seed is", seednumber))
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
  cat(proc.time() - timing)
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
