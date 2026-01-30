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
    # fit gaussian DP
    dp <- DirichletProcessGaussian(dat)
    dp <- Fit(dp, n_iter)
    
    pis    <- dp$weights                     
    mus    <- dp$clusterParameters[[1]]         
    Sigmas <- dp$clusterParameters[[2]]  
    ll = loglik_dp_v1(dat, pis, mus, Sigmas)
  }
  else{
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
dp_ll = function(dp_data, x, pax, mu, Sigma, data_scatter, data_sum){
  
}
loglik_cond_dp = function(dp_data, mu, Sigma, data_scatter, data_sum){
  
}
#---------------------- new functions ----------------------------------

ll = list()
bic = list()

child = "x5"
parents <- list(
  #c(),
  #c("x1"),
  #c("x2"),
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
pis    <- dp$weights  
K <- dp$numberClusters                
mus    <- dp$clusterParameters$mu         
Sigmas <- dp$clusterParameters$sig  

# scatter matrix
data_scatter = crossprod(dp_data)
data_sum = colSums(dp_data)


# TODO  in loop
pax = c("x1","x2")




# TODO 
#score = dp_ll(dp_data, x=child, pax, mus, Sigmas, data_scatter, data_sum)
# TODO remove
x = child
#----

vars <- c(x, pax)
pos  <- match(vars, colnames(dp_data))
child_parent_data  = dp_data[,pos]
child_parent_scatter = data_scatter[pos,pos]
child_parent_sum = data_sum[pos]

#TODO 
#for (k in 1:K) {}
k=1

pi = pis[k]
mu = mus[, , k]
child_parent_mu = mu[pos]
Sigma = Sigmas[, , k]
child_parent_Sigma = Sigma[pos,pos]

# TODO
ll = loglik_cond_dp(data=child_parent_data,
                    mu=child_parent_mu,
                    Sigma=child_parent_Sigma,
                    data_scatter=child_parent_scatter,
                    data_sum=child_parent_sum)

#TODO rmeove
data=child_parent_data
mu=child_parent_mu
Sigma=child_parent_Sigma
data_scatter=child_parent_scatter
data_sum=child_parent_sum
#-----







# define bic
n <- nrow(dat)
d <- ncol(dat)              
K <- dp$numberClusters        
p <- (K - 1) + K * d + K * (d * (d + 1) / 2)
bic <- -2 * ll + p * log(n)

list(ll = ll, bic = bic)





























for(i in 1:length(parents)){
  pax = parents[[i]]
  
  
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