#---------------------- membership prob matrix ----------------------
dp_membership_probs <- function(dp, burnin, L){
  n_iter = length(dp$alphaChain)
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
  for (l in 1:L) {
    clusterParams <- clusterParams_sample[[l]]
    pointsPerCluster <- pointsPerCluster_sample[[l]]
    numLabels <- length(pointsPerCluster)
    probs <- matrix(0, nrow = N, ncol = numLabels)

    for (i in seq_len(N)) {
      rowp <- pointsPerCluster *
        Likelihood(mdObj, y[i, , drop = FALSE], clusterParams)

      rowp[is.na(rowp)] <- 0
      if (all(rowp == 0)) {
          rowp <- rep_len(1, length(rowp))
      }

      probs[i, ] <- rowp
    }
    
    # keep only columns that are not all zero
    col_sums <- colSums(probs)
    probs <- probs[, col_sums > 0, drop = FALSE]

    probs_list[[l]] <- probs / rowSums(probs)
  }
  return(probs_list)
}
add_membershipp <- function(membershipp_list, membershipp, child, parents, vars = c(child, parents)) {
  membershipp_list[[length(membershipp_list) + 1]] <- list(
    membershipp = membershipp,
    child = child, 
    parents = parents, 
    vars = vars
  )
  membershipp_list
}
#----------------------  BiDAG ----------------------------------
usrscoreparameters <- function(initparam, 
                               usrpar = list(pctesttype = "bge",
                                             dp_iter = 100,
                                             dp_burnin = 30,
                                             dp_n_sample = 10,
                                             membershipp_list = NULL,
                                             am = 1, 
                                             aw = NULL, 
                                             T0scale = NULL,
                                             edgepf = 1,
                                             edgepmat = NULL
                                             )
                                       ) 
{
  
  if (is.null(usrpar$dp_iter)) {
    usrpar$dp_iter <- 100
  }
  if (is.null(usrpar$dp_burnin)) {
    usrpar$dp_burnin <- 30
  }
  if (is.null(usrpar$dp_n_sample)) {
    usrpar$dp_n_sample <- 10
  }
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
  initparam$dp_iter <- usrpar$dp_iter
  initparam$dp_burnin <- usrpar$dp_burnin
  initparam$dp_n_sample <- usrpar$dp_n_sample
  initparam$pf <- usrpar$edgepf
  initparam$am <- usrpar$am
  initparam$aw <- usrpar$aw
  initparam$T0scale <- usrpar$T0scale
  
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

  dp_scoreparam_list <- Filter(function(e) 
    all(param$labels[needed] %in% e$meta$vars), 
    param$dp_scoreparam_list
  )
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

## ---------------------------------------------------------------------------------------------------------------------------------------------------------
## Run MCMC
## ---------------------------------------------------------------------------------------------------------------------------------------------------------
set.searchspace <- function(data, dual, method, par = 1, alpha = 0.05, usrpar = list(pctesttype = "bge")) {
  start <- Sys.time()
  startspace <- NULL

  if(dual) {
    cor_mat <- cor(data)
    startspace <- dual_pc(cor_mat, nrow(data), alpha = alpha, skeleton = T)
  }
    
  if(method == "DP") {
    # dirichlet params
    dp_iter <- usrpar$dp_iter
    burnin <- usrpar$dp_burnin
    L <- usrpar$dp_n_sample
    dp_fits <- usrpar$dp_fits
    n = ncol(data)
    # prepare dirichlet gamma list
    Gamma_list <- list()
    for (f in seq_len(dp_fits)) {
          parents = parents,
          vars    = c(child, mb)
    }
    usrpar$membershipp_list = Gamma_list
    score <- scoreparameters("usr", data, usrpar = usrpar)
  }
  
  if(method == "bge") {
    score <- scoreparameters("bge", data, bgepar = list(am = par))
  }
  searchspace <- iterativeMCMC(scorepar = score, startspace = startspace, hardlimit = 14, 
                               verbose = F, scoreout = TRUE, alphainit = 0.01)
  time <- Sys.time() - start
  
  list(score = score, scoretable = searchspace$scoretable, DAG = searchspace$DAG, 
       maxorder = searchspace$maxorder, endspace = searchspace$endspace, time = time)
}
DP.partition.mcmc <- function(searchspace, alpha = 0.05, 
                               order = FALSE, burnin = 0.33, iterations = 600) {
  start <- Sys.time()
  Score <- searchspace$score
  
  if(order) {
    dp.fit <- orderMCMC(Score, MAP = FALSE, chainout = TRUE, alpha = alpha, 
                         startorder = searchspace$maxorder, scoretable = searchspace$scoretable,
                         startspace = searchspace$endspace, iterations = iterations, stepsave = 4)
  }
  else {
    dp.fit <- partitionMCMC(Score, alpha = alpha, startDAG = searchspace$DAG, 
                             scoretable = searchspace$scoretable, startspace = searchspace$endspace,
                             iterations = iterations, stepsave = 4)
  }
  toburn <- round(burnin * dp.fit$info$samplesteps)
  
  dp.fit$traceadd$incidence <- dp.fit$traceadd$incidence[-(1:toburn)]
  time <- Sys.time() - start + searchspace$time
  dp.fit$time <- as.numeric(time, units = "secs")
  
  return(dp.fit)
}
bge.partition.mcmc <- function(searchspace, alpha = 0.05, 
                               order = FALSE, burnin = 0.33, iterations = 600) {
  start <- Sys.time()
  BGEScore <- searchspace$score
  
  if(order) {
    bge.fit <- orderMCMC(BGEScore, MAP = FALSE, chainout = TRUE, alpha = alpha, 
                         startorder = searchspace$maxorder, scoretable = searchspace$scoretable,
                         startspace = searchspace$endspace, iterations = iterations, stepsave = 4)
  }
  else {
    bge.fit <- partitionMCMC(BGEScore, alpha = alpha, startDAG = searchspace$DAG, 
                             scoretable = searchspace$scoretable, startspace = searchspace$endspace,
                             iterations = iterations, stepsave = 4)
  }
  toburn <- round(burnin * bge.fit$info$samplesteps)
  
  bge.fit$traceadd$incidence <- bge.fit$traceadd$incidence[-(1:toburn)]
  time <- Sys.time() - start + searchspace$time
  bge.fit$time <- as.numeric(time, units = "secs")
  
  return(bge.fit)
}
## ---------------------------------------------------------------------------------------------------------------------------------------------------------
## Generate data
## 
## ---------------------------------------------------------------------------------------------------------------------------------------------------------
min_detectable_r <- function(N, q = 0, alpha = 0.05, power = 0.8) {
  T_star <- abs(qnorm(alpha / 2))
  z_power <- qnorm(power)
  tanh((T_star + z_power) / sqrt(N - q - 3))
}
pcor_from_cor <- function(R) {
  P <- -psolve(R)
  diag(P) <- 1
  P
}
make_detectable_truegraph <- function(truegraph, R, N,
                                      alpha = 0.05,
                                      power = 0.8) {
  p <- ncol(truegraph)
  q <- p - 2 # max conditioning set size for PC algorithm is p-2

  pcor <- pcor_from_cor(R)

  threshold <- min_detectable_r(
    N = N,
    q = q,
    alpha = alpha,
    power = power
  )

  detectable <- truegraph * 0
  edge_idx <- which(truegraph == 1, arr.ind = TRUE)

  for (e in seq_len(nrow(edge_idx))) {
    i <- edge_idx[e, 1]
    j <- edge_idx[e, 2]

    if (abs(pcor[i, j]) >= threshold) {
      detectable[i, j] <- 1
    }
  }
  dimnames(detectable) <- dimnames(truegraph)
  detectable
}
simulate_bimodal_student <- function(dag, n, bimodal_sep = 2, df = 3) {
  model1 <- corr(dag)
  model2 <- corr(dag)

  n1 <- sample(1:(n - 1), 1)
  n2 <- n - n1

  X1 <- rmvt(n1, sigma = model1$R, df = df)
  X2 <- rmvt(n2, sigma = model2$R, df = df)

  v <- rnorm(ncol(X2))
  v <- v / sqrt(sum(v^2))
  shift <- bimodal_sep * v
  X2 <- sweep(X2, 2, shift, "+")
  data <- standardize(rbind(X1, X2))
  if (is.null(colnames(data))) {
    colnames(data) <- paste0("v", seq_len(ncol(data)))
  }
  data  
}
simulate_bimodal <- function(dag, n, bimodal_sep=2) {
    model1 <- corr(dag)
    model2 <- corr(dag)

    n1 <- sample(1:(n - 1), 1)
    n2 <- n - n1
    
    X1 <- simulate(model1$B, model1$O, n1)
    X2 <- simulate(model2$B, model2$O, n2)
    
    v <- rnorm(ncol(X2))
    v <- v / sqrt(sum(v^2))
    shift <- bimodal_sep * v
    X2 <- sweep(X2, 2, shift, "+")
    data <- standardize(rbind(X1, X2))
    if (is.null(colnames(data))) {
      colnames(data) <- paste0("v", seq_len(ncol(data)))
    }
    data    
}
bimodal_err <- function(n, var, sep_sd=2) {
  sd0 <- sqrt(var)
  shift <- sep_sd * sd0
  comp <- sample(c(-1, 1), size=n, replace=TRUE)
  rnorm(n, mean=comp * shift, sd=sd0)
}

simulate_bimodal_one_node <- function(g, n, err=NULL, bimodal_sep=2) {
  # Randomly simulates data.
  # g = dag
  # n = sample size
  # err = additive error distribution
  
  model <- corr(g)
  B <- model$B
  O <- model$O

  # p = |variables|
  p <- ncol(B)
  
  # reorder B and O
  ord <- sofic_order(B)
  B <- B[ord, ord]
  O <- O[ord]
  
  # set default additive error as normal
  if (is.null(err)) {
    err <- function(n, var) rnorm(n, 0, sqrt(var))
  }

  # source nodes with at least one child
  source_parents_ord <- which(rowSums(B != 0) == 0 & colSums(B != 0) > 0)
  
  if (length(source_parents_ord) == 0) {
    stop("No source parent node found.")
  }
  # choose one random source parent in reordered coordinates
  chosen_ord <- sample(source_parents_ord, 1)
  
  # simulate data
  X <- matrix(0, n, p)
  for (i in 1 : p) {
    
    # linear effect
    for (j in which(B[i,] != 0)) {
      X[, i] <- X[, i] + B[i, j] * X[, j]
    }
    
    # additive error
    # additive error
    if (i == chosen_ord) {
      X[, i] <- X[, i] + bimodal_err(n, O[i], sep_sd=bimodal_sep)
    } else {
      X[, i] <- X[, i] + err(n, O[i])
    }
  }
  # reorder X
  ord <- invert_order(ord)
  X <- X[, ord]
  
  X <- standardize(X)
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("v", seq_len(ncol(X)))
  }
  return(X)  
}
