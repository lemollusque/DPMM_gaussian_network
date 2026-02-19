library(BiDAG)
library(dirichletprocess)
library(gRbase)
library(matrixStats)
library(cowplot)
library(ggplot2)

set.seed(101)

# data
N <- 100  # number of samples

x1 <- rnorm(N, mean=sample(1:10)[1], sd=1)  
x2 <- rnorm(N, mean=sample(1:10)[1], sd=1)
x3 <- rnorm(N, mean=sample(1:10)[1], sd=1)
x4 <- 1.2 * x1 - 0.8 * x2 + rnorm(N)

data <- data.frame(x1, x2, x3, x4)

#---------------------- functions ----------------------------------
dp_membership_probs <- function(dp) {
  y <- dp$data
  N <- nrow(y)
  clusterParams <- dp$clusterParameters
  numLabels <- dp$numberClusters
  mdObj <- dp$mixingDistribution
  pointsPerCluster <- dp$pointsPerCluster
  probs <- matrix(0, nrow = N, ncol = numLabels)
  for (i in seq_len(N)) {
    probs[i, 1:numLabels] <- pointsPerCluster * 
      dirichletprocess:::Likelihood.mvnormal(mdObj, 
                                             y[i,, drop = FALSE], 
                                             clusterParams)
  }
  probs <- probs / rowSums(probs)
  return(probs)
}
#----------------------  BiDAG ----------------------------------
usrscoreparameters <- function(initparam, 
                               usrpar = list(pctesttype = "bge",
                                             membershipp = NULL,
                                             am = 1, 
                                             aw = NULL, 
                                             T0scale = NULL,
                                             edgepf = 1
                                             )
                                       ) 
{
  if (is.null(usrpar$membershipp)) stop("Gamma (membershipp) is missing")
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
  initparam$pf <- usrpar$edgepf
  
  mu0 <- numeric(initparam$n)
  T0 <- diag(usrpar$T0scale, initparam$n, initparam$n)
  K = ncol(usrpar$membershipp)
  Nk <- numeric(K)
  means <- vector("list", K)
  TN <- vector("list", K)
  awpN <- numeric(K)
  constscorefact <- numeric(K)
  muN <- vector("list", K)
  SigmaN <- vector("list", K)
  for (k in  1:K){
    weightvector = usrpar$membershipp[,k]
    Nk[k] <- sum(weightvector)
    forcov <- cov.wt(initparam$data, wt = weightvector, method = "ML")
    covmatk <- forcov$cov * Nk[k]
    means[[k]] <- forcov$center
    TN[[k]] <- T0 + covmatk + 
      ((usrpar$am * Nk[k])/(usrpar$am + Nk[k])) * 
      (mu0 - means[[k]]) %*% t(mu0 - means[[k]])
    awpN[k] = usrpar$aw + Nk[k]
    constscorefact[k] =  (1/2) * log(usrpar$am/(usrpar$am + Nk[k]))
    muN[[k]] <- (Nk[k] * means[[k]] + usrpar$am * mu0)/(Nk[k] + usrpar$am)
    SigmaN[[k]] <- TN[[k]]/(awpN[k] - initparam$n - 1)
  }
  
  N <- sum(Nk)
  initparam$K <- K
  initparam$means <- means
  initparam$TN <- TN
  initparam$awpN <- awpN
  initparam$muN <- muN
  initparam$SigmaN <- SigmaN
  
  initparam$scoreconstvec <- numeric(initparam$n)
  for (j in (1:initparam$n)) {
    awp <- usrpar$aw - initparam$n + j
    initparam$scoreconstvec[j] <- -(N/2) * log(pi) + sum(constscorefact) - K*lgamma(awp/2) + 
      sum(lgamma((awp + Nk)/2)) + K*((awp + j - 1)/2) * log(usrpar$T0scale) - 
      j * log(initparam$pf)
  }
  
  initparam
}
usrDAGcorescore <- function (j, parentnodes, n, param) {
  K=param$K
  TN <- param$TN
  awpN <- param$awpN
  scoreconstvec <- param$scoreconstvec
  lp <- length(parentnodes)
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
  
  corescore
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

# perform DPMM on all parents
dp_data = scale(data) 
n_iter = 10

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
dp <- DirichletProcessMvnormal(dp_data, g0Priors)

# scoring
set.seed(12) # to test multiple clusters
dp <- Fit(dp, 10)
print(dp$numberClusters)
Gamma <- dp_membership_probs(dp)
usr_score_param <- BiDAG::scoreparameters(scoretype = "usr", 
                             data = dp_data, 
                             usrpar = list(pctesttype = "bge",
                                           membershipp = Gamma,
                                           am = alpha_mu, 
                                           aw = alpha_w, 
                                           T0scale = t,
                                           edgepf = 1
                             )
)

########################### score a DAG (check equivalence) ##################
vars  <- c("x1","x2","x3","x4")

A_12 <- matrix(c(
  0, 1, 0, 0,
  0, 0, 0, 0,
  0, 0, 0, 0,
  0, 0, 0, 0
), nrow = 4, byrow = TRUE,
dimnames = list(vars, vars))

A_21 <- matrix(c(
  0, 0, 0, 0,
  1, 0, 0, 0,
  0, 0, 0, 0,
  0, 0, 0, 0
), nrow = 4, byrow = TRUE,
dimnames = list(vars, vars))

dags <- list(
  A_12,
  A_21
)

res <- test_dag_score_equivalence(usr_score_param, dags)
res

A1 <- matrix(c(
  0,1,0,1,
  0,0,0,1,
  0,0,0,0,
  0,0,0,0
), 4, byrow=TRUE,
dimnames=list(vars, vars))

A2 <- matrix(c(
  0,1,0,1,
  0,0,0,0,
  0,0,0,0,
  0,1,0,0
), 4, byrow=TRUE,
dimnames=list(vars, vars))

A3 <- matrix(c(
  0,0,0,1,
  1,0,0,1,
  0,0,0,0,
  0,0,0,0
), 4, byrow=TRUE,
dimnames=list(vars, vars))

A4 <- matrix(c(
  0,0,0,0,
  1,0,0,1,
  0,0,0,0,
  1,0,0,0
), 4, byrow=TRUE,
dimnames=list(vars, vars))

A5 <- matrix(c(
  0,1,0,0,
  0,0,0,0,
  0,0,0,0,
  1,1,0,0
), 4, byrow=TRUE,
dimnames=list(vars, vars))


dags <- list(
  A1,
  A2,
  A3,
  A4,
  A5
)


res <- test_dag_score_equivalence(usr_score_param, dags)
res


############################### compare all dags ##############################
truegraph <- matrix(c(
  0,0,0,1,
  0,0,0,1,
  0,0,0,0,
  0,0,0,0
), 4, byrow=TRUE,
dimnames=list(vars, vars))

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
