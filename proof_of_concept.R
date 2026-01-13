setwd("C:/Users/BASTIAN/Desktop/code/git_DPMM_gaussian_network/DPMM_gaussian_network")

library(dirichletprocess)
library(gRbase)
library(BiDAG)
library(matrixStats)
library(questionr)
library(ggpubr)
library(cowplot)
library(pcalg)
library(rstan)
library(bridgesampling)

source("Fourier_fns.R")
source("BayesStanFns.R")
source("sampling_fns.R")

set.seed(101)
lambda <- 1
dual <- F  # use dualPC
n <- 4  # number of nodes
N <- 100  # number of samples

myDAG <- pcalg::randomDAG(n, prob = 0.5, lB = 1, uB = 2) 
trueDAG <- as(myDAG, "matrix")
truegraph <- 1*(trueDAG != 0)
data <- Fou_nldata(truegraph, N, lambda = lambda, noise.sd = 1, standardize = T) 

fit_joint_dp <- function(y, X, iter=200) {
  Z <- scale(cbind(X, y))              # standardize jointly is often best
  dp <- DirichletProcessMvnormal(Z)
  dp <- Fit(dp, iter)
  dp
}


get_logl = function(obs, multinomial=T, iteration=500){
  obs = scale(obs)
  if(multinomial){
    dp = DirichletProcessMvnormal(obs)
  }
  else{
    dp = DirichletProcessGaussian(obs)
  }
  dp = Fit(dp, iteration)
  
  densities = Predictive(dp$mixingDistribution, obs)
  logd = log(densities)
  logl = sum(logd)
  
  return(logl)
}


cond_loglik_from_joint <- function(dp, y, X) {
  #Standardize joint data 
  Z  <- scale(cbind(X, y))
  d  <- ncol(as.matrix(X))
  Xs <- Z[, 1:d, drop = FALSE]
  ys <- Z[, d + 1]
  
  N <- nrow(Xs)
  
  #Extract fitted mixture parameters 
  pis    <- dp$weights                      # length K
  mus    <- dp$clusterParameters$mu         # list, each length d+1
  Sigmas <- dp$clusterParameters$sig        # list, each (d+1)x(d+1)
  
  K <- length(pis)
  stopifnot(length(mus)/ncol(Z) == K, length(Sigmas)/(ncol(Z)^2) == K)
  
  logp <- numeric(N)
  
  #Loop over observations 
  for (i in 1:N) {
    
    xi <- matrix(Xs[i, ], nrow = d)
    
    log_w_unnorm <- numeric(K)
    log_comp_y   <- numeric(K)
    
    # Loop over mixture components 
    for (k in 1:K) {
      
      mu  <- mus[,,k]
      Sig <- Sigmas[,,k]
      
      ## Partition mean and covariance
      mu_x <- matrix(mu[1:d], nrow = d)
      mu_y <- mu[d + 1]
      
      S_xx <- Sig[1:d,     1:d,     drop = FALSE]
      S_xy <- Sig[1:d,     d + 1,   drop = FALSE]
      S_yx <- Sig[d + 1,   1:d,     drop = FALSE]
      S_yy <- Sig[d + 1,   d + 1]
      
      #x-dependent mixture weights 
      log_px <- mvtnorm::dmvnorm(
        x     = t(xi),
        mean  = as.numeric(mu_x),
        sigma = S_xx,
        log   = TRUE
      )
      
      log_w_unnorm[k] <- log(pis[k]) + log_px
      
      # Conditional Gaussian y | x 
      invSxx <- solve(S_xx)
      
      m <- mu_y + as.numeric(S_yx %*% invSxx %*% (xi - mu_x))
      V <- as.numeric(S_yy - S_yx %*% invSxx %*% S_xy)
      
      ## numerical safety
      V <- max(V, 1e-10)
      
      log_comp_y[k] <- dnorm(
        ys[i],
        mean = m,
        sd   = sqrt(V),
        log  = TRUE
      )
    }
    
    #Mixture log-density ------------------------------------------
    log_w <- log_w_unnorm - matrixStats::logSumExp(log_w_unnorm)
    
    logp[i] <- matrixStats::logSumExp(log_w + log_comp_y)
  }
  
  #Return conditional log score -----------------------------------
  sum(logp)
}




##############################################################
# List all DAGs with n nodes
all.dags <- list()
adj <- matrix(0, nrow = n, ncol = n)
dag.counter <- 0
all.comb <- rep(list(c(0,1)), n*(n-1))
all.comb <- expand.grid(all.comb)  # all combinations outside of diagonal of adjacency matrix

for(i in 1:nrow(all.comb)) {
  adj[col(adj)!=row(adj)] <- as.numeric(all.comb[i, ])
  
  if(is.DAG(adj)) {
    dag.counter <- dag.counter + 1
    all.dags[[dag.counter]] <- adj
  }
}

# Compute true posterior (and save all scores)
true.post <- rep(NA, dag.counter)
parent.scores <- data.frame(sets = character(0), newscore = character(0))  
parent.scores <- rep(list(parent.scores), n)

for(k in 1:dag.counter) {
  print(k)
  dag <- all.dags[[k]]
  curr_score <- 0
  
  for(x in 1:n) {
    print(x)
    set <- parent.scores[[x]]$sets
    pax <- dag[,x]
    pax.str <- paste(pax, collapse = "")  # parents of x
    check <- is.na(match(set, pax.str))  # check if parent set is already there
    
    if(all(check)) {  # new parent set
      
      if(sum(pax) == 0) {  # Compute local score
        print("no pax")
        dp <- DirichletProcessGaussian(scale(data[,x]))
        dp <- Fit(dp, 200)
        loc_score <- sum(log(Predictive(dp$mixingDistribution, scale(data[,x]))))
      }
      else {
        print("pax")
        dp <- DirichletProcessMvnormal(scale(cbind(data[ ,which(pax==1)], data[,x])))
        dp <- Fit(dp, 200)
        loc_score <- cond_loglik_from_joint(dp, data[,x], data[ ,which(pax==1)])
      }
      parent.scores[[x]][length(set)+1, ] <- c(pax.str, loc_score)
    }
    
    else {  # fetch score from parent.scores
      ind <- which(check == F)  # index of pax in set
      loc_score <- as.numeric(parent.scores[[x]]$newscore[ind])
    }
    curr_score <- curr_score + loc_score  # build score
  }
  true.post[k] <- curr_score
}

true.p <- exp(true.post - logSumExp(true.post))


#compute shd for every graph
shd <- rep(NA, dag.counter)
count_edges <- rep(NA, dag.counter)
for(k in 1:dag.counter) {
  print(k)
  dag <- all.dags[[k]]
  camparison = compareDAGs(dag, truegraph)
  shd[k] = camparison["SHD"]
  count_edges[k] = sum(dag)
}
shd.order =  order(shd, decreasing = F) 
a <- true.post[shd.order]
b <- true.p[shd.order]

plot(a, type="l")
plot(b, type="l")


