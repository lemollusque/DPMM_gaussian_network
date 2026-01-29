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
library(ggplot2)

set.seed(101)
N <- 100  # number of samples

x1 <- rnorm(N, mean=sample(1:10)[1], sd=1)  
x2 <- rnorm(N, mean=sample(1:10)[1], sd=1)
x3 <- rnorm(N, mean=sample(1:10)[1], sd=1)
x4 <- rnorm(N, mean=sample(1:10)[1], sd=1)
x5 <- 1.2 * x1 - 0.8 * x2 + rnorm(N)

data <- data.frame(x1, x2, x3, x4, x5)

loglik_cond_component_dp <- function(dp_data, mu, Sigma) { 
  y <- dp_data[,1] #child node is in the first column
  X <- as.matrix(dp_data[,-1, drop=FALSE]) #parents
  
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
  
  #entry (i,k) = log pi_k + log p_k(y_i | x_i)
  ll_mat <- matrix(NA_real_, nrow = n, ncol = K)
  for (k in 1:K) {
    ll_mat[, k] <- log(pis[k]) + loglik_cond_component_dp(dp_data, mus[, , k], Sigmas[, , k])
  }
  # total log-likelihood
  m <- apply(ll_mat, 1, max)
  sum(m + log(rowSums(exp(ll_mat - m))))
}

loglik_dp <- function(dp_data, pis, mus, Sigmas) {
  n <- nrow(dp_data)
  K <- length(pis)
  
  #entry (i,k) = log pi_k + log p_k(y_i | x_i)
  ll_mat <- matrix(NA_real_, nrow = n, ncol = K)
  for (k in 1:K) {
    ll_mat[, k] <- log(pis[k]) + dnorm(dp_data, mean = mus[, , k], sd = sqrt(Sigmas[, , k]), log = TRUE)
  }
  # total log-likelihood
  m <- apply(ll_mat, 1, max)
  sum(m + log(rowSums(exp(ll_mat - m))))
}

dp_ll = function(dat, n_iter){
  #first column = child node
  #remaining columns = parents
  if (ncol(dat)<2){
    # fit gaussian DP
    dp <- DirichletProcessGaussian(dat)
    dp <- Fit(dp, n_iter)
    
    pis    <- dp$weights                     
    mus    <- dp$clusterParameters[[1]]         
    Sigmas <- dp$clusterParameters[[2]]  
    ll = loglik_dp(dat, pis, mus, Sigmas)
  }
  else{
    # Fit multivariate normal
    dp <- DirichletProcessMvnormal(dat)
    dp <- Fit(dp, n_iter)
    
    pis    <- dp$weights                     
    mus    <- dp$clusterParameters$mu         
    Sigmas <- dp$clusterParameters$sig   
    ll = loglik_cond_dp(dat, pis, mus, Sigmas)
  }
  # define bic
  n <- nrow(dat)
  d <- ncol(dat)              
  K <- dp$numberClusters        
  p <- (K - 1) + K * d + K * (d * (d + 1) / 2)
  bic <- -2 * ll + p * log(n)
  
  list(ll = ll, bic = bic)
  
}



ll = list()
bic = list()

parents <- list(
  c(),
  c("x1"),
  c("x2"),
  c("x1","x2"),
  c("x1","x2","x3"),
  c("x1","x2","x4"),
  c("x1","x2","x3","x4")
)
for(i in 1:length(parents)){
  set.seed(101)
  pax = parents[[i]]
  dp_data = scale(data[ ,c("x5" , pax)]) 
  
  
  # loglikelihood
  score = dp_ll(dp_data, 20)
  ll[[i]] <- score$ll
  bic[[i]] <- score$bic
  
}

labs <- sapply(parents, function(p) if (length(p)==0) "none" else paste(p, collapse=","))

par(mfrow = c(3,1))

plot(1:length(parents), ll, pch=19, xaxt="n",
     xlab="Parent set", ylab="LL")
axis(1, seq_len(length(parents)), labels=labs, las=2)

plot(1:length(parents), bic, pch=19, xaxt="n",
     xlab="Parent set", ylab="BIC")
axis(1, seq_len(length(parents)), labels=labs, las=2)

post <- -0.5 * unlist(bic)
p <- exp(post - logSumExp(post))
plot(1:length(parents), p, pch=19, xaxt="n",
     xlab="Parent set", ylab="softmax")
axis(1, seq_len(length(parents)), labels=labs, las=2)