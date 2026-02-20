#---------------------- membership prob matrix ----------------------
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
test_compare_dp_vs_bge <- function(usr_score_param,
                                   bge_score_param,
                                   dags,tol = 1e-4,
                                   verbose = TRUE) {
  if (is.null(names(dags)) || any(names(dags) == "")) {
    names(dags) <- paste0("DAG_", seq_along(dags))
  }
  
  # score all dags
  scores <- unlist(lapply(dags, function(A) BiDAG::DAGscore(usr_score_param, A)))
  bge_scores <- unlist(lapply(dags, function(A) BiDAG::DAGscore(bge_score_param, A)))
  
  diffs = scores-bge_scores
  equal <- max(abs(diffs)) < tol
  
  if (verbose) {
    cat("Total dp scores:\n")
    print(scores)
    cat("Total bge scores:\n")
    print(bge_scores)
    cat("Equal?", equal, "\n\n")
  }
  equal
}