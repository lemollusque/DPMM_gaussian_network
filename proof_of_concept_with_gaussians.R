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


#---------------------- functions ----------------------------------
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

dp_ll = function(dp, child, pax){ 
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
      
      ll_mat[, k] <- log(pis[k]) + loglik_cond_component_dp(child_parent_data, child_parent_mu, child_parent_Sigma)
      
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
  list_bic = list()
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
        list_bic[[idx]] <- score$bic
        idx <- idx + 1
      } else {
        cmb <- combn(others, k, simplify = FALSE)
        for (s in cmb) {
          score = dp_ll(dp, child, c(pax, s))
          list_bic[[idx]] <- score$bic
          idx <- idx + 1
        }
      }
    }
  }
  bics = unlist(list_bic)
  log_evid <- -0.5 * bics
  
  m <- max(log_evid)
  log_mean_evid <- m + log(mean(exp(log_evid - m)))
  
  -2 * log_mean_evid
  
}
}
#----------------------  functions ----------------------------------

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
dp_list <- list()
vars  <- c("x1","x2","x3","x4","x5")
for (child in vars){
  parents <- vars[vars != child]
  dp_data = scale(data[,c(child, parents)]) 
  n_iter = 20
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