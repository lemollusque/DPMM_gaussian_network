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

x1 <- rnorm(N, mean=sample(1:10)[1], sd=sample(1:10)[1]) # TODO check wiith 0 ,1 
x2 <- rnorm(N, mean=sample(1:10)[1], sd=sample(1:10)[1])
x3 <- rnorm(N, mean=sample(1:10)[1], sd=sample(1:10)[1])
x4 <-rnorm(N, mean=sample(1:10)[1], sd=sample(1:10)[1])
x5 <- 1.2 * x1 - 0.8 * x2 + rnorm(N)

data <- data.frame(x1, x2, x3, x4, x5)
# 
# trueDAG <- matrix(0, nrow = ncol(data), ncol = ncol(data),
#                   dimnames = list(names(data), names(data)))
# # parents 
# trueDAG["x1", "x5"] <- 1
# trueDAG["x2", "x5"] <- 1
# 
# # add DAGs with more parents
# all.dags <- list()
# all.dags[[1]] <- trueDAG
# 
# d = trueDAG
# d["x3", "x5"] <- 1
# all.dags[[2]] <- d
# 
# d = trueDAG
# d["x4", "x5"] <- 1
# all.dags[[3]] <- d
# 
# d = trueDAG
# d["x3", "x5"] <- 1
# d["x4", "x5"] <- 1
# all.dags[[4]] <- d




loglik_cond_component_dp <- function(dp_data, mu, Sigma) { 
  
  y <- dp_data[,1]
  X <- as.matrix(dp_data[,-1, drop=FALSE])
  
  mu_y <- mu[1]
  mu_x <- mu[-1]
  
  S_yy <- Sigma[1,1]
  S_yx <- Sigma[1,-1, drop=FALSE]  
  S_xy <- Sigma[-1,1, drop=FALSE]   
  S_xx <- Sigma[-1,-1, drop=FALSE] 
  
  
  # Cholesky of S_xx 
  U <- chol(S_xx)
  v <- backsolve(U, forwardsolve(t(U), t(S_yx)))
  beta <- as.numeric(v)     
  
  w <- backsolve(U, forwardsolve(t(U), S_xy))     
  quad <- as.numeric(S_yx %*% w)                 
  sigma2 <- S_yy - quad
  
  Xm <- sweep(X, 2, mu_x, "-")    
  mean_cond <- mu_y + as.numeric(Xm %*% beta)
  
  resid <- y - mean_cond
  
  ll <- -0.5 * (log(2*pi) + log(sigma2) + (resid^2)/sigma2)
  ll


}


logsumexp_vec <- function(a) {
  m <- max(a)
  m + log(sum(exp(a - m)))
}



loglik_cond_dp <- function(dp_data, pis, mus, Sigmas) {
  n <- nrow(dp_data)
  K <- length(pis)
  logpis <- log(pis)
  
  # compute n x C matrix of log(pi_c) + log p_c(y|x)
  ll_mat <- matrix(NA_real_, n, K)
  for (k in 1:K) {
    ll_mat[,k] <- logpis[k] + loglik_cond_component_dp(dp_data, mus[,,k], Sigmas[,,k])
  }
  
  # row-wise logsumexp
  ll_i <- apply(ll_mat, 1, logsumexp_vec)
  sum(ll_i)
}

ll = list()

parents = list()
parents[[1]] = c("x1", "x2")
parents[[2]] = c("x1", "x2", "x3")
parents[[3]] = c("x1", "x2", "x4")
parents[[4]] = c("x1", "x2", "x3", "x4")
for(i in 1:length(parents)){
  pax = parents[[i]]
  dp_data = scale(data[ ,c("x5" , pax)]) 
  dp <- DirichletProcessMvnormal(dp_data)
  dp <- Fit(dp, 1000)
  
  pis    <- dp$weights                     
  mus    <- dp$clusterParameters$mu         
  Sigmas <- dp$clusterParameters$sig   
  
  
  ll[[i]] <- loglik_cond_dp(dp_data, pis, mus, Sigmas)
}
  
plot(1:length(parents), ll)


# compare means and covs
pax = parents[[4]]
dp_data = scale(data[ ,c("x5" , pax)])


dp <- DirichletProcessMvnormal(dp_data)
dp <- Fit(dp, 200)
                   
mus    <- dp$clusterParameters$mu         
Sigmas <- dp$clusterParameters$sig 

Delta <- Sigmas[,,1] - cov(dp_data)

cov(dp_data)
mus
Sigmas
Delta































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




##############################################################
false.post = true.post[2:dag.counter]
false.p =   exp(false.post - logSumExp(false.post))
#compute shd for every graph
shd <- rep(NA, dag.counter)
count_edges <- rep(NA, dag.counter)
for(k in 2:dag.counter) {
  print(k)
  dag <- all.dags[[k]]
  camparison = compareDAGs(dag, truegraph)
  shd[k] = camparison["SHD"]
  count_edges[k] = sum(dag)
}
shd.order =  order(shd, decreasing = F) 
a <- false.post[shd.order]
b <- false.p[shd.order]

plot(a, type="l")
plot(b, type="l")
