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

# data
N <- 100  # number of samples

x1 <- rnorm(N, mean=sample(1:10)[1], sd=1)  
x2 <- rnorm(N, mean=sample(1:10)[1], sd=1)
x3 <- rnorm(N, mean=sample(1:10)[1], sd=1)
x4 <- rnorm(N, mean=sample(1:10)[1], sd=1)
x5 <- 1.2 * x1 - 0.8 * x2 + rnorm(N)

data <- data.frame(x1, x2, x3, x4, x5)


#---------------------- old functions ----------------------------------
loglik_cond_component_dp_v1 <- function(dp_data, mu, Sigma) { 
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


loglik_cond_dp_v1 <- function(dp_data, pis, mus, Sigmas) {
  n <- nrow(dp_data)
  K <- length(pis)
  
  #entry (i,k) = log pi_k + log p_k(y_i | x_i)
  ll_mat <- matrix(NA_real_, nrow = n, ncol = K)
  for (k in 1:K) {
    ll_mat[, k] <- log(pis[k]) + loglik_cond_component_dp_v1(dp_data, mus[, , k], Sigmas[, , k])
  }
  # total log-likelihood
  m <- apply(ll_mat, 1, max)
  sum(m + log(rowSums(exp(ll_mat - m))))
}

loglik_dp_v1 <- function(dp_data, pis, mus, Sigmas) {
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

dp_ll_v1 = function(dat, n_iter){
  #first column = child node
  #remaining columns = parents
  if (ncol(dat)<2){
    set.seed(101)
    # fit gaussian DP
    dp <- DirichletProcessGaussian(dat)
    dp <- Fit(dp, n_iter)
    
    pis    <- dp$weights                     
    mus    <- dp$clusterParameters[[1]]         
    Sigmas <- dp$clusterParameters[[2]]  
    ll = loglik_dp_v1(dat, pis, mus, Sigmas)
  }
  else{
    set.seed(101)
    # Fit multivariate normal
    dp <- DirichletProcessMvnormal(dat)
    dp <- Fit(dp, n_iter)
    
    pis    <- dp$weights                     
    mus    <- dp$clusterParameters$mu         
    Sigmas <- dp$clusterParameters$sig   
    ll = loglik_cond_dp_v1(dat, pis, mus, Sigmas)
  }
  # define bic
  n <- nrow(dat)
  d <- ncol(dat)              
  K <- dp$numberClusters        
  p <- (K - 1) + K * d + K * (d * (d + 1) / 2)
  bic <- -2 * ll + p * log(n)
  
  list(ll = ll, bic = bic)
  
}
#---------------------- end old functions ----------------------------------


#---------------------- new functions ----------------------------------
loglik_cond_component_dp = function(n, mu, Sigma, data_scatter, data_sum){
  mu_y <- mu[1]
  mu_x <- mu[-1]
  
  S_yy <- Sigma[1,1]
  S_yx <- Sigma[1,-1, drop=FALSE]   
  S_xx <- Sigma[-1,-1, drop=FALSE] 
  
  # quadratic part of the trace
  C = data_scatter - data_sum %*% t(mu) - mu %*% t(data_sum) + n* mu %*% t(mu)
  C_yy <- C[1,1]
  C_yx <- C[1,-1, drop=FALSE]    
  C_xx <- C[-1,-1, drop=FALSE] 
  
  # v = Syy - Syx Sxx^-1 Sxy
  L = t(base::chol(S_xx))
  w <- forwardsolve(L, t(S_yx))
  u <- backsolve(t(L), w)
  v = as.numeric(S_yy - S_yx %*% u)
  
  # sscss = Syx Sxx^-1 Cxx Sxx^-1 Sxy
  css <- C_xx %*% u  
  w <- forwardsolve(L, css)
  u <- backsolve(t(L), w)
  sscss = as.numeric(S_yx %*% u)
  
  # ssc = Syx Sxx^−1 Cxy
  w <- forwardsolve(L, t(C_yx))
  u <- backsolve(t(L), w)
  ssc = as.numeric(S_yx %*% u)
  
  # conditional quadratic term q
  q = (C_yy - 2*ssc + sscss)/v
  
  -(n/2)*log(2*base::pi*v) - (1/2)*q
}


loglik_component_dp <- function(n, mu, Sigma, data_scatter, data_sum){
  C <- data_scatter - 2*mu*data_sum + n*mu^2
  -(n/2)*log(2*base::pi*Sigma) - 0.5*C/Sigma
}


dp_ll = function(dp, child, pax){ 
  dp_data =dp$data
  pis <- dp$weights  
  K <- length(pis)                
  mus    <- dp$clusterParameters$mu         
  Sigmas <- dp$clusterParameters$sig  
  
  
  
  n <- nrow(dp_data)
  vars <- c(child, pax)
  pos  <- match(vars, colnames(dp_data))
  
  # initiate ll
  ll = 0
  
  for (k in 1:K) {
    # cluster params
    mu_k = mus[, , k]
    child_parent_mu = mu_k[pos]
    Sigma_k = Sigmas[, , k]
    child_parent_Sigma = Sigma_k[pos,pos, drop=FALSE]
    
    # data in cluster k
    data_k = dp_data[dp$clusterLabels == k,, drop=FALSE]
    n_k <- nrow(data_k)
    
    ## data scatter 
    data_scatter = crossprod(data_k)
    data_sum = colSums(data_k)
    child_parent_scatter = data_scatter[pos,pos, drop=FALSE]
    child_parent_sum = data_sum[pos]
    
    if (length(pax) < 1){
      # compute complete-data ll for cluster k
      ll_component = loglik_component_dp(n=n_k,
                                         mu=child_parent_mu,
                                         Sigma=child_parent_Sigma,
                                         data_scatter=child_parent_scatter,
                                         data_sum=child_parent_sum)
      ll <- ll + n_k * log(pis[k]) + ll_component
    }
    else{
      # compute conditional complete-data ll for cluster k
      ll_cond_component = loglik_cond_component_dp(n=n_k,
                                                   mu=child_parent_mu,
                                                   Sigma=child_parent_Sigma,
                                                   data_scatter=child_parent_scatter,
                                                   data_sum=child_parent_sum)
      ll <- ll + n_k * log(pis[k]) + ll_cond_component
    }
  }
  
  # define bic
  d <- length(vars)             
  p <- (K - 1) + K * d + K * (d * (d + 1) / 2)
  bic <- -2 * ll + p * log(n)
  
  list(ll = ll, bic = bic)
}
#---------------------- new functions ----------------------------------
child = "x5"
parents <- list(
  c(),
  c("x1"),
  c("x2"),
  c("x1","x2"),
  c("x1","x2","x3"),
  c("x1","x2","x4"),
  c("x1","x2","x3","x4")
)
possible_parents= sort(unique(unlist(parents)))
# perform DPMM on all parents
dp_data = scale(data[,c(child, possible_parents)]) 
n_iter = 200
dp <- DirichletProcessMvnormal(dp_data)
dp <- Fit(dp, n_iter)

ll = list()
bic = list()
ll_old = list()
bic_old = list()
for(i in 1:length(parents)){
  pax = parents[[i]]
  
  # loglikelihood
  score = dp_ll(dp, child, pax)
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