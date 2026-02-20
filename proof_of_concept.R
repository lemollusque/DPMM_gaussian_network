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
dp_membership_probs <- function(dp, n_iter, burnin, L){
  # define samples to take
  post_idx <- seq(burnin + 1, n_iter)
  thin_idx <- round(seq(length(post_idx), 1, length.out = L))
  sample_idx <- post_idx[thin_idx]
  
  y <- dp$data
  N <- dp$n
  mdObj <- dp$mixingDistribution # const if updatePrior = FALSE
  clusterParams_sample = dp$clusterParametersChain[sample_idx]
  weightsChain_sample = dp$weightsChain[sample_idx]
  pointsPerCluster_sample <- lapply(weightsChain_sample, function(w) w * N) 
  
  
  probs_list <- vector("list", length(clusterParams_sample))
  for (l in  1:L){
    clusterParams = clusterParams_sample[[l]]
    pointsPerCluster <- pointsPerCluster_sample[[l]]
    numLabels <- length(pointsPerCluster)
    probs <- matrix(0, nrow = N, ncol = numLabels)
    for (i in seq_len(N)) {
      probs[i, 1:numLabels] <- pointsPerCluster * 
        dirichletprocess:::Likelihood.mvnormal(mdObj, 
                                               y[i,, drop = FALSE], 
                                               clusterParams)
    }
    probs_list[[l]] <- probs / rowSums(probs)
  }
  return(probs_list)
}
add_membershipp <- function(membershipp_list, membershipp, child, parents) {
  membershipp_list[[length(membershipp_list) + 1]] <- list(
    membershipp = membershipp,
    child = child, 
    parents = parents, 
    vars = c(child, parents) 
  )
  membershipp_list
}
#----------------------  BiDAG ----------------------------------
usrscoreparameters <- function(initparam, 
                               usrpar = list(pctesttype = "bge",
                                             membershipp_list = NULL,
                                             am = 1, 
                                             aw = NULL, 
                                             T0scale = NULL,
                                             edgepf = 1,
                                             edgepmat = NULL
                               )
) 
{
  if (is.null(usrpar$membershipp_list) || length(usrpar$membershipp_list) == 0)
    stop("membershipp_list is missing or empty")
  if (is.null(usrpar$edgepf)) {
    usrpar$edgepf <- 1
  }
  if (is.null(usrpar$am)) {
    usrpar$am <- 1
  }
  if (is.null(usrpar$aw)) {
    usrpar$aw <- initparam$n + usrpar$am + 1
  }
  if (is.null(usrpar$T0scale)) {
    usrpar$T0scale <- usrpar$am * (usrpar$aw - initparam$n - 1)/(usrpar$am + 1)
  }
  if (is.null(usrpar$edgepmat)) {
    initparam$logedgepmat <- NULL
  }
  else {
    initparam$logedgepmat <- log(usrpar$edgepmat)
  }
  
  initparam$pf <- usrpar$edgepf
  initparam$am <- usrpar$am
  initparam$aw <- usrpar$aw
  
  # only depending on n
  mu0 <- numeric(initparam$n)
  T0 <- diag(usrpar$T0scale, initparam$n, initparam$n)
  
  # loop over all DPs
  dp_membershipp_list<- usrpar$membershipp_list
  n_dp <- length(dp_membershipp_list)
  initparam$dp_scoreparam_list <- vector("list", n_dp)
  for (d in 1:n_dp){
    membershipp_list = dp_membershipp_list[[d]]$membershipp
    L <- length(membershipp_list)
    scoreparam_list <- vector("list", L)
    for (l in 1:L) {
      membershipp = membershipp_list[[l]]
      K = ncol(membershipp)
      Nk <- numeric(K)
      means <- vector("list", K)
      TN <- vector("list", K)
      awpN <- numeric(K)
      constscorefact <- numeric(K)
      for (k in  1:K){
        weightvector = membershipp[,k]
        Nk[k] <- sum(weightvector)
        forcov <- cov.wt(initparam$data, wt = weightvector, method = "ML")
        covmatk <- forcov$cov * Nk[k]
        means[[k]] <- forcov$center
        TN[[k]] <- T0 + covmatk + 
          ((usrpar$am * Nk[k])/(usrpar$am + Nk[k])) * 
          (mu0 - means[[k]]) %*% t(mu0 - means[[k]])
        awpN[k] = usrpar$aw + Nk[k]
        constscorefact[k] =  (1/2) * log(usrpar$am/(usrpar$am + Nk[k]))
      }
      
      N <- sum(Nk)
      
      scoreconstvec <- numeric(initparam$n)
      for (j in (1:initparam$n)) {
        awp <- usrpar$aw - initparam$n + j
        scoreconstvec[j] <- -(N/2) * log(pi) + sum(constscorefact) - K*lgamma(awp/2) + 
          sum(lgamma((awp + Nk)/2)) + K*((awp + j - 1)/2) * log(usrpar$T0scale) - 
          j * log(initparam$pf)
      }
      
      # save score params for DP_list[l]
      scoreparam_list[[l]] <- list(
        K = K,
        TN = TN,
        awpN = awpN,
        scoreconstvec = scoreconstvec
      )
    }
    initparam$dp_scoreparam_list[[d]] <- list(
      meta = list(
        child   = dp_membershipp_list[[d]]$child,
        parents = dp_membershipp_list[[d]]$parents,
        vars    = dp_membershipp_list[[d]]$vars
      ),
      scores = scoreparam_list
    )
  }
  initparam
}
usrDAGcorescore <- function (j, parentnodes, n, param) {
  # not depending on cluster param
  lp <- length(parentnodes)
  
  # extract needed score parameters
  needed <- c(j, parentnodes)
  dp_scoreparam_list = Filter(function(e) all(param$labels[needed] %in% e$meta$vars), 
                              param$dp_scoreparam_list)
  # put all score parameters needed in one list.
  scoreparam_list = unlist(lapply(dp_scoreparam_list, `[[`, "scores"), recursive = FALSE)
  n_score = length(scoreparam_list)
  corescore_list  <- numeric(n_score)
  for (s in 1:n_score){
    scoreparam = scoreparam_list[[s]]
    K=scoreparam$K
    TN <- scoreparam$TN
    awpN <- scoreparam$awpN
    scoreconstvec <- scoreparam$scoreconstvec
    awpNd2 <- (awpN - n + lp + 1)/2
    A <- sapply(TN, function(m) m[j, j])
    
    switch(as.character(lp), 
           `0` = {
             corescore <- scoreconstvec[lp + 1] - sum(awpNd2 * log(A))
           }, 
           `1` = {
             D <- sapply(TN, function(m) m[parentnodes, parentnodes])
             logdetD <- log(D)
             B <- sapply(TN, function(m) m[j, parentnodes])
             logdetpart2 <- log(A - B^2/D)
             corescore <- scoreconstvec[lp + 1] - sum(awpNd2 * logdetpart2) - 
               sum(logdetD)/2
             if (!is.null(param$logedgepmat)) {
               corescore <- corescore - param$logedgepmat[parentnodes, 
                                                          j]
             }
           }, 
           `2` = {
             D <- lapply(TN, function(m) m[parentnodes, parentnodes, drop = FALSE])
             detD <- sapply(D, function(m) BiDAG:::dettwobytwo(m))
             logdetD <- log(detD)
             B <- lapply(TN, function(m) m[j, parentnodes, drop = FALSE])
             logdetpart2 <- vapply(seq_along(D), function(k) {
               Ak <- A[k]
               Bk <- matrix(B[[k]], nrow = 1)            
               Mk <- D[[k]] - (t(Bk) %*% Bk) / Ak        
               log(BiDAG:::dettwobytwo(Mk)) + log(Ak) - logdetD[k]
             }, numeric(1))
             corescore <- scoreconstvec[lp + 1] - sum(awpNd2 * logdetpart2) - 
               sum(logdetD)/2
             if (!is.null(param$logedgepmat)) {
               corescore <- corescore - sum(param$logedgepmat[parentnodes, 
                                                              j])
             }
           }, 
           {
             D <- lapply(TN, function(m) as.matrix(m[parentnodes, parentnodes, drop = FALSE]))
             choltemp <- lapply(D, function(m) chol(m))
             logdetD <- vapply(seq_along(TN), function(k) {
               2 * sum(log(diag(choltemp[[k]])))
             }, numeric(1))
             B <- lapply(TN, function(m) m[j, parentnodes, drop = FALSE])
             logdetpart2 <- vapply(seq_along(TN), function(k) {
               val <- log(A[k] - sum(backsolve(choltemp[[k]], t(B[[k]]), transpose=TRUE)^2))
             }, numeric(1))
             corescore <- scoreconstvec[lp + 1] - sum(awpNd2 * logdetpart2) - 
               sum(logdetD)/2
             if (!is.null(param$logedgepmat)) {
               corescore <- corescore - sum(param$logedgepmat[parentnodes, 
                                                              j])
             }
           })
    corescore_list[s] = corescore
  }
  mean(corescore_list)
}
#----------------------  test functions ----------------------------------
test_dag_score_equivalence <- function(usr_score_param,
                                       dags, 
                                       tol = 1e-4,
                                       verbose = TRUE) {
  if (is.null(names(dags)) || any(names(dags) == "")) {
    names(dags) <- paste0("DAG_", seq_along(dags))
  }
  
  # score all dags
  scores <- unlist(lapply(dags, function(A) BiDAG::DAGscore(usr_score_param, A)))
  
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
# replace BIDAG functions
unlockBinding("usrscoreparameters", asNamespace("BiDAG"))
assign("usrscoreparameters", usrscoreparameters, envir = asNamespace("BiDAG"))
lockBinding("usrscoreparameters", asNamespace("BiDAG"))

unlockBinding("usrDAGcorescore", asNamespace("BiDAG"))
assign("usrDAGcorescore", usrDAGcorescore, envir = asNamespace("BiDAG"))
lockBinding("usrDAGcorescore", asNamespace("BiDAG"))


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
n_iter = 100
burnin = 30
L = 10 # sample to take

Gamma_list <- list()
vars  <- c("x1","x2","x3","x4")
for (child in vars){
  parents <- vars[vars != child]
  dp_data = scaled_data[,c(child, parents)]
  dp <-  DirichletProcessMvnormal(dp_data, g0Priors)
  dp <- Fit(dp, n_iter)
  
  Gamma_sample <- dp_membership_probs(dp, n_iter, burnin, L)
  Gamma_list <- add_membershipp(Gamma_list, Gamma_sample, child=child, parents=parents)
}

# scoring
usr_score_param <- BiDAG::scoreparameters(scoretype = "usr", 
                                          data = scaled_data, 
                                          usrpar = list(pctesttype = "bge",
                                                        membershipp_list = Gamma_list,
                                                        am = alpha_mu, 
                                                        aw = alpha_w, 
                                                        T0scale = t,
                                                        edgepf = 1
                                          )
)
############################### compare all dags ##############################
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

# Compute true posterior (and save all scores)
true.post <- rep(NA, dag.counter)

for(d in 1:dag.counter) {
  print(d)
  dag <- all.dags[[d]]
  true.post[d] = BiDAG::DAGscore(usr_score_param, dag)
}
# Order true posterior
true.order <- order(true.post, decreasing = T)
true.post <- true.post[true.order]
true.p <- exp(true.post - logSumExp(true.post))
all.dags <- all.dags[true.order]


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


par(mfrow = c(2,1))
plot(a, type="l", ylab="posterior")
plot(b, type="l", ylab="softmax")




top5_idx <- order(true.post, decreasing = T)[1:5]
all.dags[top5_idx]
true.post[top5_idx]
