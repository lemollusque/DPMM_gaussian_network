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



#----------------------  overwrite functions ----------------------------------
# replace BIDAG functions
unlockBinding("usrscoreparameters", asNamespace("BiDAG"))
assign("usrscoreparameters", usrscoreparameters, envir = asNamespace("BiDAG"))
lockBinding("usrscoreparameters", asNamespace("BiDAG"))

unlockBinding("usrDAGcorescore", asNamespace("BiDAG"))
assign("usrDAGcorescore", usrDAGcorescore, envir = asNamespace("BiDAG"))
lockBinding("usrDAGcorescore", asNamespace("BiDAG"))
#----------------------  end overwrite functions ------------------------------
## ---------------------------------------------------------------------------------------------------------------------------------------------------------
## Load and prepare the data
## ---------------------------------------------------------------------------------------------------------------------------------------------------------
dataname = "gaussian"
set.seed(12)

# data
N <- 100  # number of samples

x1 <- rnorm(N, mean=sample(1:10)[1], sd=1)  
x2 <- rnorm(N, mean=sample(1:10)[1], sd=1)
x3 <- rnorm(N, mean=sample(1:10)[1], sd=1)
x4 <- 1.2 * x1 - 0.8 * x2 + rnorm(N)

data <- data.frame(x1, x2, x3, x4)

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
n_iter = 50
burnin = 30
L = 10 # sample to take

cormat <- cor(scaled_data)
pc.skel = pcalg::pc(suffStat = list(C = cormat, 
                          n = N), indepTest = pcalg::gaussCItest, alpha = 0.05, 
          labels = colnames(scaled_data), skel.method = "stable", 
          verbose = FALSE)

g <- pc.skel@graph
startspace <- 1 * (graph2m(g))


Gamma_list <- list()
vars  <- c("x1","x2","x3","x4")
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

# score
usr_score_param <- BiDAG::scoreparameters(scoretype = "usr", 
                                          data = scaled_data, 
                                          usrpar = list(pctesttype = "bge",
                                                        dp_iter = n_iter,
                                                        dp_burnin = burnin,
                                                        dp_n_sample = L,
                                                        membershipp_list = Gamma_list,
                                                        am = alpha_mu, 
                                                        aw = alpha_w, 
                                                        T0scale = t,
                                                        edgepf = 1
                                          )
)

# add parents
new_space = startspace
new_space[3,4] = 1
new_space[2,1] = 1
new_space[1,2] = 1

update_score_param <- function(param, space){
  for (i in 1:length(param$dp_scoreparam_list)){
    # deactivate all
    param$dp_scoreparam_list[[i]]$meta$active = F
  }
  # initiate new probability matrices
  Gamma_list <- list()
  for (child in colnames(space)){
    parents <- names(which(space[ , child] == 1))
    idx <- which(vapply(param$dp_scoreparam_list, function(e) {
      identical(e$meta$child, child) &&
        setequal(e$meta$parents, parents)
    }, logical(1)))
    
    # activate
    if (length(idx) > 0) {
      if (length(idx) > 1) stop("multiple DP with same child parents")
      for (j in idx) {
        param$dp_scoreparam_list[[j]]$meta$active <- TRUE
      }
    }
    
    # create new DP
    if (length(idx) == 0) {
      dp_data = param$data[,c(child, parents)]
      if (length(parents) == 0){
        dp <-  DirichletProcessGaussian(dp_data,
                                        alphaPriors = c(1, 20))
      }
      else{
        n_col = ncol(dp_data)
        g0Priors <- list(
          mu0    = rep(0, n_col),
          Lambda = diag(n_col) / param$T0scale,   # T = (1/t) I
          kappa0 = param$am,
          nu     = param$aw
        )
        
        dp <-  DirichletProcessMvnormal(dp_data, g0Priors)
      }
      dp <- Fit(dp, param$dp_iter)
      
      Gamma_sample <- dp_membership_probs(dp, param$dp_burnin, param$dp_n_sample)
      # add meta data in list
      Gamma_list <- add_membershipp(Gamma_list, 
                                    Gamma_sample, 
                                    child=child, 
                                    parents=parents, 
                                    active=TRUE)
      
    }
    
  }
  new_score = BiDAG::scoreparameters(scoretype = "usr", 
                                     data = param$data, 
                                     usrpar = list(pctesttype = "bge",
                                                   membershipp_list = Gamma_list,
                                                   am = param$am, 
                                                   aw = param$aw, 
                                                   T0scale = param$T0scale,
                                                   edgepf = 1
                                     )
  )
  # add new computed params to existing param
  param$dp_scoreparam_list <- append(param$dp_scoreparam_list, new_score$dp_scoreparam_list)
  param
}

new_score_param = update_score_param(usr_score_param, new_space)

########################### score a DAG (check equivalence) ##################
A_12 <- matrix(c(
  0, 1, 0, 0,
  0, 0, 0, 0,
  0, 0, 0, 0,
  0, 0, 0, 0
), nrow = 4, byrow = TRUE,
dimnames = list(vars, vars))

A_21 <- matrix(c(
  0, 0, 0, 0,
  1, 0, 0, 0,
  0, 0, 0, 0,
  0, 0, 0, 0
), nrow = 4, byrow = TRUE,
dimnames = list(vars, vars))

dags <- list(
  A_12,
  A_21
)

# only true if both  dags are in search space
test_dag_score_equivalence(new_score_param, dags)


