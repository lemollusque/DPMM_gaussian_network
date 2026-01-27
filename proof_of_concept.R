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


loglik_cond_component_dp <- function(dp_data, mu, Sigma) { 
  
  y <- dp_data[,1] #child node is in the first column
  X <- as.matrix(dp_data[,-1, drop=FALSE])
  
  mu_y <- mu[1]
  mu_x <- mu[-1]
  
  S_yy <- Sigma[1,1]
  S_yx <- Sigma[1,-1, drop=FALSE]  
  S_xy <- Sigma[-1,1, drop=FALSE]   
  S_xx <- Sigma[-1,-1, drop=FALSE] 
  
  
  # Cholesky of S_xx 
  L = t(base::chol(S_xx))
  
  #conditional distr (mu and Sigma)
  w <- forwardsolve(L, t(X) - mu_x)
  z = backsolve(t(L), w)
  mu_cond = drop(mu_y + S_yx %*% z)
  
  W <- forwardsolve(L, S_xy)
  Z = backsolve(t(L), W)
  S_cond = S_yy - as.numeric(S_yx %*% Z)
  
  # log density
  dnorm(y, mean = mu_cond, sd = sqrt(S_cond), log = TRUE)
  
}


loglik_cond_dp <- function(dp_data, pis, mus, Sigmas) {
  n <- nrow(dp_data)
  K <- length(pis)
  
  logpis <- log(pis)
  
  #entry (i,k) = log pi_k + log p_k(y_i | x_i)
  ll_mat <- matrix(NA_real_, nrow = n, ncol = K)
  for (k in 1:K) {
    ll_mat[, k] <- logpis[k] + loglik_cond_component_dp(dp_data, mus[, , k], Sigmas[, , k])
  }
  
  # total log-likelihood: sum_i logsumexp_k ll_mat[i,k]
  m <- apply(ll_mat, 1, max)
  sum(m + log(rowSums(exp(ll_mat - m))))
}

loglik_dp <- function(dp_data, pis, mus, Sigmas) {
  n <- nrow(dp_data)
  K <- length(pis)
  
  logpis <- log(pis)
  
  #entry (i,k) = log pi_k + log p_k(y_i | x_i)
  ll_mat <- matrix(NA_real_, nrow = n, ncol = K)
  for (k in 1:K) {
    ll_mat[, k] <- logpis[k] + dnorm(dp_data, mean = mus[, , k], sd = sqrt(Sigmas[, , k]), log = TRUE)
  }
  
  # total log-likelihood: sum_i logsumexp_k ll_mat[i,k]
  m <- apply(ll_mat, 1, max)
  sum(m + log(rowSums(exp(ll_mat - m))))
}

dp_ll = function(dat, n_iter){
  #first column = child node
  # remaining parents
  if (ncol(dat)<2){
    
    dp <- DirichletProcessGaussian(dat)
    dp <- Fit(dp, n_iter)
    
    pis    <- dp$weights                     
    mus    <- dp$clusterParameters[[1]]         
    Sigmas <- dp$clusterParameters[[2]]  
    ll = loglik_dp(dat, pis, mus, Sigmas)
    
    
    
  }
  else{
    dp <- DirichletProcessMvnormal(dat)
    dp <- Fit(dp, n_iter)
    
    pis    <- dp$weights                     
    mus    <- dp$clusterParameters$mu         
    Sigmas <- dp$clusterParameters$sig   
    ll = loglik_cond_dp(dat, pis, mus, Sigmas)
  }
  
  n <- nrow(dat)
  d <- ncol(dat)              
  K <- length(dp$numberClusters)         
  p <- (K - 1) + K * d + K * (d * (d + 1) / 2)
  bic <- -2 * ll + p * log(n)
  
  
  list(ll = ll, bic = bic)
  
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
    set <- parent.scores[[x]]$sets
    pax <- dag[,x]
    pax.str <- paste(pax, collapse = "")  # parents of x
    check <- is.na(match(set, pax.str))  # check if parent set is already there
    
    if(all(check)) {  # new parent set
      dp_data = scale(data[ ,c(x , pax)])
      
      ll = dp_ll(dp_data, 1000)
      loc_score <- ll$bic
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