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

dp_ll_v1 = function(dp, child, pax){ 
  dp_data=dp$data
  pis <- dp$weights  
  K <- length(pis)                
  mus    <- dp$clusterParameters$mu         
  Sigmas <- dp$clusterParameters$sig  
  
  
  
  n <- nrow(dp_data)
  vars <- c(child, pax)
  d <- length(vars)    
  pos  <- match(vars, colnames(dp_data))
  child_parent_data = dp_data[,pos]
  
  # initiate ll
  ll_mat <- matrix(NA_real_, nrow = n, ncol = K)
  for (k in 1:K) {
    # cluster params
    mu_k = mus[, , k]
    child_parent_mu = mu_k[pos]
    Sigma_k = Sigmas[, , k]
    child_parent_Sigma = Sigma_k[pos,pos, drop=FALSE]
    
    if (length(pax) < 1){
      ll_mat[, k] <- log(pis[k]) + dnorm(child_parent_data, mean = child_parent_mu, sd = sqrt(child_parent_Sigma), log = TRUE)
    }
    else{
      
      ll_mat[, k] <- log(pis[k]) + loglik_cond_component_dp_v1(child_parent_data, child_parent_mu, child_parent_Sigma)
      
    }
    
    # total log-likelihood
    m <- apply(ll_mat, 1, max)
    ll = sum(m + log(rowSums(exp(ll_mat - m))))
  }
  # define bic         
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
    data_k = dp_data
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
  
  # define BIC
  d <- length(vars)             
  p <- (K - 1) + K * d + K * (d * (d + 1) / 2)
  bic <- -2 * ll + p * log(n)
  
  list(ll = ll, bic = bic)
}

#---------------------- end new functions ----------------------------------
child = "x5"
parents_list <- list(
  c(),
  c("x1"),
  c("x2"),
  c("x1","x2"),
  c("x1","x2","x3"),
  c("x1","x2","x4"),
  c("x1","x2","x3","x4")
)
possible_parents= sort(unique(unlist(parents_list)))
# perform DPMM on all parents
dp_data = scale(data[,c(child, possible_parents)]) 
n_iter = 200
dp <- DirichletProcessMvnormal(dp_data)
dp <- Fit(dp, n_iter)

ll = list()
bic = list()
ll_old = list()
bic_old = list()
for(i in 1:length(parents_list)){
  pax = parents_list[[i]]
  
  # loglikelihood
  score = dp_ll_v1(dp, child, pax)
  ll_old[[i]] <- score$ll
  bic_old[[i]] <- score$bic
  
  
  # loglikelihood
  score = dp_ll(dp, child, pax)
  ll[[i]] <- score$ll
  bic[[i]] <- score$bic
}

# plot results
labs <- sapply(parents_list, function(p) if (length(p)==0) "none" else paste(p, collapse=","))

par(mfrow = c(3,1))

plot(1:length(parents_list), ll, pch=19, xaxt="n",
     xlab="Parent set", ylab="LL")
axis(1, seq_len(length(parents_list)), labels=labs, las=2)

plot(1:length(parents_list), bic, pch=19, xaxt="n",
     xlab="Parent set", ylab="BIC")
axis(1, seq_len(length(parents_list)), labels=labs, las=2)

post <- -0.5 * unlist(bic)
p <- exp(post - logSumExp(post))
plot(1:length(parents_list), p, pch=19, xaxt="n",
     xlab="Parent set", ylab="softmax")
axis(1, seq_len(length(parents_list)), labels=labs, las=2)


########################### fit dps and then average ##################

#---------------------------- avg dp functions ----------------------------
add_dp <- function(dp_list, dp, child, parents) {
  dp_list[[length(dp_list) + 1]] <- list(
    dp = dp,
    child = child,
    parents = parents,
    vars = c(child, parents) 
  )
  dp_list
}

find_dps <- function(dp_list, child, parents) {
  needed <- c(child, parents)
  Filter(function(e) all(needed %in% e$vars), dp_list)
}

average_dp_ll = function(dp_list, child, pax){ 
  avg_ll = list()
  idx = 1
  for (i in 1:length(dp_list)){
    dp = dp_list[[i]]$dp
    
    vars = dp_list[[i]]$vars
    required <- unique(c(child, pax))
    others <- setdiff(vars, required)
    
    # choose any subset of the remaining vars
    for (k in 0:length(others)) {
      if (k == 0) {
        # loglikelihood
        score = dp_ll(dp, child, pax)
        avg_ll[[idx]] <- score$bic
        idx <- idx + 1
      } else {
        cmb <- combn(others, k, simplify = FALSE)
        for (s in cmb) {
          score = dp_ll(dp, child, c(pax, s))
          avg_ll[[idx]] <- score$bic
          idx <- idx + 1
        }
      }
    }
  }
  mean(unlist(avg_ll))
}

#---------------------------- end avg dp functions ----------------------------
dp_list <- list()
vars  <- c("x1","x2","x3","x4","x5")
for (child in vars){
  parents <- vars[vars != child]
  dp_data = scale(data[,c(child, parents)]) 
  n_iter = 200
  dp <- DirichletProcessMvnormal(dp_data)
  dp <- Fit(dp, n_iter)
  dp_list <- add_dp(dp_list, dp, child=child, parents=parents)
}

parents_list <- list(
  c(),
  c("x1"),
  c("x2"),
  c("x1","x2"),
  c("x1","x2","x3"),
  c("x1","x2","x4"),
  c("x1","x2","x3","x4")
)

avg_score <- list()
for(i in 1:length(parents_list)){
  pax = parents_list[[i]]
  
  hits <- find_dps(dp_list, child, pax)
  
  # bic
  avg_score[[i]] = average_dp_ll(hits, child, pax)
}

# plot avg BIC
labs <- sapply(parents_list, function(p) if (length(p)==0) "none" else paste(p, collapse=","))
par(mfrow = c(2,1))

plot(1:length(parents_list), avg_score, pch=19, xaxt="n",
     xlab="Parent set", ylab="AVG BIC")
axis(1, seq_len(length(parents_list)), labels=labs, las=2)
axis(1, seq_len(length(parents_list)), labels=labs, las=2)

post <- -0.5 * unlist(avg_score)
p <- exp(post - logSumExp(post))
plot(1:length(parents_list), p, pch=19, xaxt="n",
     xlab="Parent set", ylab="softmax")
axis(1, seq_len(length(parents_list)), labels=labs, las=2)