library(dirichletprocess)
library(gRbase)
library(BiDAG)
library(matrixStats)
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

vars  <- c("x1","x2","x3","x4")
colnames(data) = vars


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
    
  }
  
  # total log-likelihood
  m <- apply(ll_mat, 1, max)
  ll = sum(m + log(rowSums(exp(ll_mat - m))))
  
  # define bic         
  k_cond <- length(pax)*K
  bic <- -2 * ll + k_cond * log(n)
  
  list(ll = ll, bic = bic)
}

cov_ll = function(dp_data, child, pax){   
  n <- nrow(dp_data)
  vars <- c(child, pax)
  d <- length(vars)    
  pos  <- match(vars, colnames(dp_data))
  
  mu    <- colMeans(dp_data)        
  Sigma <- cov(dp_data) * (n - 1) / n
  
  child_parent_data = dp_data[,pos]
  child_parent_mu = mu[pos]
  child_parent_Sigma = Sigma[pos,pos, drop=FALSE]
  
  # initiate ll
  ll_mat <- matrix(NA_real_, nrow = n, ncol = 1)
  
  if (length(pax) < 1){
    ll_mat[, 1] <- dnorm(child_parent_data, mean = child_parent_mu, sd = sqrt(child_parent_Sigma), log = TRUE)
  }
  else{
    ll_mat[, 1] <- loglik_cond_component_dp(child_parent_data, child_parent_mu, child_parent_Sigma)
  }
  
  # total log-likelihood
  m <- apply(ll_mat, 1, max)
  ll = sum(m + log(rowSums(exp(ll_mat - m))))
  
  # define bic         
  k_cond <- length(pax)
  bic <- -2 * ll + k_cond * log(n)
  
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
  for (i in 1:length(dp_list)){
    dp = dp_list[[i]]$dp
    
    # loglikelihood
    score = dp_ll(dp, child, pax)
    list_bic[[i]] <- score$bic
  }
  bics = unlist(list_bic)
  log_evid <- -0.5 * bics
  
  m <- max(log_evid)
  log_mean_evid <- m + log(mean(exp(log_evid - m)))
  
  -2 * log_mean_evid
  
}

score_a_dag = function(adj_mat, dp_list){
  vars = colnames(adj_mat)
  local_scores <- numeric(length(vars))
  names(local_scores) <- vars
  for (child in vars){
    pax = names(adj_mat[,child][which(adj_mat[,child] == 1)])
    
    hits <- find_dps(dp_list, child, pax)
    
    # bic
    local_scores[child] = average_dp_ll(hits, child, pax)
  }
  sum(local_scores)
}
#----------------------  test functions ----------------------------------
test_dag_score_equivalence <- function(dags, dp_list,
                                       tol = 1e-4,
                                       verbose = TRUE) {
  if (is.null(names(dags)) || any(names(dags) == "")) {
    names(dags) <- paste0("DAG_", seq_along(dags))
  }
  
  # score all dags
  scores <- unlist(lapply(dags, function(A) score_a_dag(A, dp_list)))
  
  # pairwise differences
  diffs <- outer(scores, scores, FUN = "-")
  
  # pass/fail: all equal within tol?
  all_equal <- max(abs(diffs)) < tol
  
  if (verbose) {
    cat("Total scores:\n")
    print(scores)
    cat("\nMax |pairwise diff|:", max(abs(diffs)), "\n")
    cat("All equal?", all_equal, "\n\n")
  }
  
  all_equal
}
#----------------------  functions ----------------------------------

# List all DAGs with n nodes

dp_list <- list()
for (child in vars){
  parents <- vars[vars != child]
  dp_data = scale(data[,c(child, parents)]) 
  n_iter = 200
  dp <- DirichletProcessMvnormal(dp_data)
  dp <- Fit(dp, n_iter)
  dp_list <- add_dp(dp_list, dp, child=child, parents=parents)
}
# List all DAGs with n nodes
all.dags <- list()
adj <- matrix(0, nrow = n, ncol = n,
              dimnames=list(vars, vars))
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


scores <- rep(NA, dag.counter)
for(d in 1:dag.counter) {
  print(d)
  dag <- all.dags[[d]]
  scores[d] = score_a_dag(dag, dp_list)
}


true.p <- exp(-0.5*scores - logSumExp(-0.5*scores))


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
a <- scores[shd.order]
b <- true.p[shd.order]


par(mfrow = c(2,1))
plot(a, type="l", ylab="score")
plot(b, type="l", ylab="softmax")
