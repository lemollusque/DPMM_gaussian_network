#---------------------- membership prob matrix ----------------------
dp_membership_probs <- function(fit, L = 500) {
  Y <- as.matrix(fit$data)
  N <- nrow(Y)
  p <- ncol(Y)
  M <- length(fit$mean)

  sample_idx <- unique(round(seq(1, M, length.out = min(L, M))))
  probs_list <- vector("list", length(sample_idx))

  for (ll in seq_along(sample_idx)) {
    m <- sample_idx[ll]
    mu <- as.matrix(fit$mean[[m]])
    Sigma <- fit$sigma2[[m]]
    w <- as.numeric(fit$probs[[m]])

    K <- length(w)
    log_probs <- matrix(-Inf, N, K)

    for (k in seq_len(K)) {
      if (p == 1) {
        sig2 <- as.numeric(Sigma)[k]
        log_probs[, k] <-
          log(w[k]) +
          dnorm(
            Y[, 1],
            mean = mu[k],
            sd = sqrt(sig2),
            log = TRUE
          )

      } else {
        Sk <- Sigma[, , k, drop = FALSE][, , 1]
        log_probs[, k] <-
          log(w[k]) +
          mvtnorm::dmvnorm(
            Y,
            mean = mu[k, ],
            sigma = Sk,
            log = TRUE
          )
      }
    }


    row_max <- apply(log_probs, 1, max)
    probs <- exp(log_probs - row_max)
    probs[!is.finite(probs)] <- 0

    zero_rows <- rowSums(probs) == 0
    if (any(zero_rows)) {
      probs[zero_rows, ] <- 1 / K
    }

    keep <- colSums(probs) > 0
    probs <- probs[, keep, drop = FALSE]

    probs_list[[ll]] <- probs / rowSums(probs)
  }
  probs_list
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
# helper function
BGeaugment <- function(sigma, mu, N, n, am, aw, logedgepmat, pf) {
  mu0 <- numeric(length(mu))
  T0scale <- am * (aw - n - 1)/(am + 1)
  T0 <- diag(T0scale, length(mu), length(mu))
  TN <- T0 + sigma + ((am * N)/(am + N)) * (mu0 - mu) %*% t(mu0 - mu)
  awpN <- aw + N
  constscorefact <- -(N/2) * log(pi) + (1/2) * log(am/(am + N))
  muN <- (N * mu + am * mu0)/(N + am)
  SigmaN <- TN/(awpN - n - 1)
  scoreconstvec <- numeric(length(mu))
  for (j in (1:length(mu))) {
    awp <- aw - n + j
    scoreconstvec[j] <- constscorefact - lgamma(awp/2) + 
      lgamma((awp + N)/2) + ((awp + j - 1)/2) * log(T0scale) - j * log(pf)
  }
  localparam <- list()
  localparam$type <- "bge"
  localparam$TN <- TN
  localparam$awpN <- awpN
  localparam$n <- n
  localparam$scoreconstvec <- scoreconstvec
  localparam$DBN <- FALSE
  localparam$MDAG <- FALSE
  localparam$logedgepmat <- logedgepmat
  return(localparam)
}
usrscoreparameters <- function(initparam, 
                               usrpar = list(pctesttype = "bge",
                                             dp_mcmc = list(niter = 5000, nburn = 3000, model="LS"),
                                             dp_prior = list(strength = 1, discount = 0),
                                             dp_n_sample = 10,
                                             dp_fits = 1,
                                             membershipp_list = NULL,
                                             am = 1, 
                                             aw = NULL, 
                                             T0scale = NULL,
                                             edgepf = 1,
                                             edgepmat = NULL
                                             )
                                       ) 
{
  if (is.null(usrpar$dp_prior)) {
    usrpar$dp_prior <- list(strength = 1, discount = 0, model="LS")
  }
  if (is.null(usrpar$dp_mcmc)) {
    usrpar$dp_mcmc <- list(niter = 5000, nburn = 3000)
  }
  if (is.null(usrpar$dp_fits)) {
    usrpar$dp_fits <- 1
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
  initparam$dp_prior <- usrpar$dp_prior
  initparam$dp_mcmc <- usrpar$dp_mcmc
  initparam$dp_fits <- usrpar$dp_fits
  initparam$dp_n_sample <- usrpar$dp_n_sample
  initparam$pf <- usrpar$edgepf
  initparam$am <- usrpar$am
  initparam$aw <- usrpar$aw
  initparam$T0scale <- usrpar$T0scale
    
  # only depending on n
  n = initparam$n
  mu0 <- numeric(n)
  T0 <- diag(usrpar$T0scale, n, n)
  
  # loop over all DPs
  dp_membershipp_list<- usrpar$membershipp_list
  n_dp <- length(dp_membershipp_list)
  initparam$dp_scoreparam_list <- vector("list", n_dp)
  mu_sum <- numeric(n)
  cov_sum <- matrix(0, n, n)
  for (d in 1:n_dp){
    membershipp_list = dp_membershipp_list[[d]]$membershipp
    L <- length(membershipp_list)
    scoreparam_list <- vector("list", L)
    for (l in 1:L) {
      membershipp = membershipp_list[[l]]
      K = ncol(membershipp)
      Nk <- numeric(K)
      means <- vector("list", K)
      covs <- vector("list", K)
      TN <- vector("list", K)
      awpN <- numeric(K)
      constscorefact <- numeric(K)
      for (k in  1:K){
        weightvector = membershipp[,k]
        Nk[k] <- sum(weightvector)
        forcov <- cov.wt(initparam$data, wt = weightvector, method = "ML")
        covs[[k]] <- forcov$cov
        covmatk <- forcov$cov * Nk[k]
        means[[k]] <- forcov$center
        TN[[k]] <- T0 + covmatk + 
          ((usrpar$am * Nk[k])/(usrpar$am + Nk[k])) * 
          (mu0 - means[[k]]) %*% t(mu0 - means[[k]])
        awpN[k] = usrpar$aw + Nk[k]
        constscorefact[k] =  (1/2) * log(usrpar$am/(usrpar$am + Nk[k]))
      }
      
      N <- sum(Nk)
      clusterWeights <- Nk/N

      scoreconstvec <- numeric(n)
      for (j in (1:n)) {
        awp <- usrpar$aw - n + j
        scoreconstvec[j] <- -(N/2) * log(pi) + sum(constscorefact) - K*lgamma(awp/2) + 
          sum(lgamma((awp + Nk)/2)) + K*((awp + j - 1)/2) * log(usrpar$T0scale) 
      }

      # for bge params, we need the average mean and cov across clusters, weighted by cluster weights
      mu_bar <- Reduce("+", Map(function(mu, w) w * mu, means, clusterWeights))
      cov_bar <- Reduce("+", Map(function(S, w) w * S, covs, clusterWeights))
      mu_sum <- mu_sum + mu_bar
      cov_sum <- cov_sum + cov_bar

      # save score params for DP_list[l]
      scoreparam_list[[l]] <- list(
        K = K,
        TN = TN,
        awpN = awpN,
        scoreconstvec = scoreconstvec,
        clusterWeights = clusterWeights
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

  # set up bge params with averaged means and covs across DPs
  N <- nrow(initparam$data)
  means  <- mu_sum / (n_dp * usrpar$dp_n_sample)
  covmat  <- cov_sum / (n_dp * usrpar$dp_n_sample)
  covmat  <- covmat * N
  bgeinitparam <- BGeaugment(covmat, means, N, n, usrpar$am, usrpar$aw, initparam$logedgepmat, initparam$pf)
  initparam$bgeinitparam <- bgeinitparam

  initparam
}
usrDAGcorescore <- function (j, parentnodes, n, param) {
  BiDAG:::DAGcorescore(j, parentnodes, n, param$bgeinitparam)
}
targetCoreScore <- function (j, parentnodes, n, param, l) {
  # not depending on cluster param
  lp <- length(parentnodes)
  # extract needed score parameters
  needed <- c(j, parentnodes)
  
  scoreparam_list = unlist(lapply(param$dp_scoreparam_list, `[[`, "scores"), recursive = FALSE)
  scoreparam <- scoreparam_list[[l]]
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
  corescore
}
DPscoreDAG <- function(param, dag) {
  n <- ncol(dag)
  scoreparam_list = unlist(lapply(param$dp_scoreparam_list, `[[`, "scores"), recursive = FALSE)
  L <- length(scoreparam_list)
  dag_scores <- numeric(L)

  for(l in 1:L){
    curr_score <- 0
    for(x in 1:n) {
      parents <- which(dag[, x] == 1)
      loc_score <- targetCoreScore(x, parents, n, param, l)    
      curr_score <- curr_score + loc_score  # build score
    }
    dag_scores[l] <- curr_score
  }
  logMeanExpLogs (dag_scores)
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
  scores <- unlist(lapply(dags, function(A) DPscoreDAG(usr_score_param, A)))
  
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
  scores <- unlist(lapply(dags, function(A) DPscoreDAG(usr_score_param, A)))
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
set.searchspace <- function(data, method, par = 1, alpha = 0.05, usrpar = list(pctesttype = "bge")) {
  start <- Sys.time()
  cor_mat <- cor(data)
  startspace <- dual_pc(cor_mat, nrow(data), alpha = alpha, skeleton = T)

  if(method == "DP") {
    # dirichlet params
    prior <- usrpar$dp_prior
    mcmc <- usrpar$dp_mcmc
    L <- usrpar$dp_n_sample
    dp_fits <- usrpar$dp_fits

    # compute fitspace
    if(usrpar$dp_fitspace == "pc") {
      cor_mat <- cor(data)
      pc.skel <- pcalg::skeleton(suffStat = list(C = cor_mat, 
                    n = nrow(data)), indepTest = pcalg::gaussCItest, alpha = alpha, 
                    labels = colnames(data), method = "stable", 
                    verbose = FALSE)
      g <- pc.skel@graph
      fitspace <- 1 * (graph2m(g))
    }
    if(usrpar$dp_fitspace == "dual") {
      fitspace <- startspace
    }  
    if(usrpar$dp_fitspace == "full") {
      fitspace <- 1 - diag(ncol(data))
      dimnames(fitspace) <- list(colnames(data), colnames(data))
    }

    # prepare dirichlet gamma list
    Gamma_list <- list()
    for (f in seq_len(dp_fits)) {
      for (child in colnames(fitspace)){
        parents <- names(which(fitspace[ , child] == 1))
        children <- names(which(fitspace[child, ] == 1))
        neighbour <- setdiff(unique(c(parents, children)), child)
        
        # create DP
        dp_data = data[,c(child, neighbour)]
        output <- list(out_param = TRUE, out_type = "FULL")  
        dp <- PYdensity(y = dp_data, mcmc = mcmc, prior = prior, output = output)
        
        # add meta data in list
        Gamma_sample <- dp_membership_probs(dp, L)
        Gamma_list <- add_membershipp(
          Gamma_list, 
          Gamma_sample, 
          child   = child, 
          parents = parents,
          vars    = c(child, neighbour)
        )
        
      }
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
                               order = FALSE, burnin = 0.33, iterations = 1000) {
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
  inter <- Sys.time()

  toburn <- round(burnin * dp.fit$info$samplesteps)
  dp.fit$traceadd$incidence <- dp.fit$traceadd$incidence[-(1:toburn)]
  dp.fit$trace <- dp.fit$trace[-(1:toburn)]
  ndags <- length(dp.fit$trace)

  # compute weights for all sampled unique DAGs
  dag_key <- function(dag) {
    paste(c(dim(dag), as.integer( 1 * as(dag, "matrix"))), collapse = "_")
  }
  dags <- dp.fit$traceadd$incidence
  keys <- vapply(dags, dag_key, character(1))
  unique_keys <- unique(keys)
  first_idx <- match(unique_keys, keys)
  unique_scores <- numeric(length(unique_keys))
  names(unique_scores) <- unique_keys

  for (u in seq_along(unique_keys)) {
    dag <- dags[[first_idx[u]]]

    unique_scores[u] <- DPscoreDAG(
      param = Score,
      dag = dag
    )
  }

  target_scores <- unname(unique_scores[keys])
  weights <- target_scores - dp.fit$trace

  dp.fit$weights <- weights - logSumExp(weights)
  end <- Sys.time()
  time2 <- end - inter
  time <- end - start + searchspace$time
  dp.fit$time <- as.numeric(time, units = "secs")
  dp.fit$time2 <- as.numeric(time2, units = "secs")
  return(dp.fit)
}
bge.partition.mcmc <- function(searchspace, alpha = 0.05, 
                               order = FALSE, burnin = 0.33, iterations = 1000) {
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
B_to_R <- function(B, O) {
  p <- ncol(B)
  IB <- solve(diag(p) - B)
  S <- IB %*% diag(O) %*% t(IB)
  cov2cor(S)
}
simulate_bimodal_student <- function(dag, n, bimodal_sep = 2, df = 3,
                             return_model = FALSE) {
  model1 <- corr(dag)
  model2 <- corr(dag)

  n1 <- sample(floor(0.2 * n):ceiling(0.8 * n), 1)
  n2 <- n - n1

  det1 <- make_detectable_truegraph(truegraph = dag, R = model1$R,  N = n1)
  det2 <- make_detectable_truegraph(truegraph = dag, R = model2$R,  N = n2)

  detectable <- 1 * ((det1 + det2) > 0)

  model1$B <- model1$B * detectable
  model2$B <- model2$B * detectable

  model1$R <- B_to_R(model1$B, model1$O)
  model2$R <- B_to_R(model2$B, model2$O)

  X1 <- rmvt(n1, sigma = model1$R, df = df)
  X2 <- rmvt(n2, sigma = model2$R, df = df)

  v <- rnorm(ncol(X2))
  v <- v / sqrt(sum(v^2))
  shift <- bimodal_sep * v
  X2 <- sweep(X2, 2, shift, "+")

  raw_mean = colMeans(rbind(X1, X2))
  raw_sd = apply(rbind(X1, X2), 2, sd)

  data <- standardize(rbind(X1, X2))
  if (is.null(colnames(data))) {
    colnames(data) <- paste0("v", seq_len(ncol(data)))
  }

  if (return_model) {
    return(list(
      data = data,
      model1 = model1,
      model2 = model2,
      n1 = n1,
      n2 = n2,
      raw_mean = raw_mean,
      raw_sd = raw_sd,
      detectable_truegraph = as(detectable, "matrix")
    ))
  }

  data  
}
simulate_bimodal <- function(dag, n, bimodal_sep = 2,
                             return_model = FALSE) {
  model1 <- corr(dag)
  model2 <- corr(dag)

  n1 <- sample(floor(0.2 * n):ceiling(0.8 * n), 1)
  n2 <- n - n1

  det1 <- make_detectable_truegraph(truegraph = dag, R = model1$R,  N = n1)
  det2 <- make_detectable_truegraph(truegraph = dag, R = model2$R,  N = n2)

  detectable <- 1 * ((det1 + det2) > 0)

  model1$B <- model1$B * detectable
  model2$B <- model2$B * detectable

  model1$R <- B_to_R(model1$B, model1$O)
  model2$R <- B_to_R(model2$B, model2$O)

  X1 <- simulate(model1$B, model1$O, n1)
  X2 <- simulate(model2$B, model2$O, n2)

  v <- rnorm(ncol(X2))
  v <- v / sqrt(sum(v^2))
  shift <- bimodal_sep * v
  X2 <- sweep(X2, 2, shift, "+")

  raw_mean = colMeans(rbind(X1, X2))
  raw_sd = apply(rbind(X1, X2), 2, sd)

  data <- standardize(rbind(X1, X2))
  if (is.null(colnames(data))) {
    colnames(data) <- paste0("v", seq_len(ncol(data)))
  }

  if (return_model) {
    return(list(
      data = data,
      model1 = model1,
      model2 = model2,
      n1 = n1,
      n2 = n2,
      raw_mean = raw_mean,
      raw_sd = raw_sd,
      detectable_truegraph = as(detectable, "matrix")
    ))
  }

  data
}
bimodal_err <- function(n1, n2, var, sep_sd=2) {
  sd0 <- sqrt(var)
  shift <- sep_sd * sd0
  c(
    rnorm(n1, mean = -shift, sd = sd0),
    rnorm(n2, mean =  shift, sd = sd0)
  )
}

simulate_bimodal_one_node <- function(g, n, err=NULL, bimodal_sep=2,
                             return_model = FALSE) {
  # Randomly simulates data.
  # g = dag
  # n = sample size
  # err = additive error distribution
  
  model <- corr(g)

  n1 <- sample(floor(0.2 * n):ceiling(0.8 * n), 1)
  n2 <- n - n1

  det <- make_detectable_truegraph(truegraph = g, R = model$R,  N = n)
  detectable <- 1 * (det > 0)
  model$B <- model$B * detectable
  model$R <- B_to_R(model$B, model$O)

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
    if (i == chosen_ord) {
      X[, i] <- X[, i] + bimodal_err(n1, n2, O[i], sep_sd=bimodal_sep)
    } else {
      X[, i] <- X[, i] + err(n, O[i])
    }
  }
  # save node
  chosen_original <- ord[chosen_ord]
  bimodal_node_name = colnames(X)[chosen_original]

  # reorder X
  ord <- invert_order(ord)
  X <- X[, ord]
  
  raw_mean = colMeans(X)
  raw_sd = apply(X, 2, sd)

  X <- standardize(X)
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("v", seq_len(ncol(X)))
  }

  if (return_model) {
    return(list(
      data = X,
      model = model,
      n1 = n1,
      n2 = n2,
      raw_mean = raw_mean,
      raw_sd = raw_sd,
      bimodal_node = chosen_original,
      bimodal_node_name = bimodal_node_name,
      detectable_truegraph = as(detectable, "matrix")
      ))
  }

  return(X)  
}
## ---------------------------------------------------------------------------------------------------------------------------------------------------------
## Non interventional effects
## Bestie equivalent functions
## ---------------------------------------------------------------------------------------------------------------------------------------------------------
DP_sample_node_beta <- function(j, parentNodes, dataParams, sample = TRUE) {
  # chose the DP score parameters corresponding to the child node j
  child_name <- colnames(dataParams$data)[j]
  d_idx <- which(vapply(
    dataParams$dp_scoreparam_list,
    function(x) x$meta$child == child_name,
    logical(1)
  ))

  # Chose one DP draw for the scoring:
  # If several DP fits exist for the same child, choose one
  d <- sample(d_idx, 1)
  scoreparam_list <- dataParams$dp_scoreparam_list[[d]]$scores
  # Choose one posterior DP membership draw
  if (sample) {
    l <- sample(seq_along(scoreparam_list), 1)
  } else {
    l <- 1
  }
  scoreparam <- scoreparam_list[[l]]

  # find clusters and cluster weights
  K <- scoreparam$K
  beta_by_cluster <- vector("list", K)

  # same as Bestie:::DAGparametersCore but for each cluster
  for (k in seq_len(K)) {
    TNk <- scoreparam$TN[[k]]
    df <- scoreparam$awpN[k] - dataParams$n + length(parentNodes) + 1
    R11 <- TNk[parentNodes, parentNodes, drop = FALSE]
    R12 <- TNk[parentNodes, j, drop = FALSE]
    R11inv <- solve(R11)
    mb <- R11inv %*% R12
    divisor <- TNk[j, j] - t(R12) %*% mb
    sigma <- as.numeric(divisor / df) * R11inv
    if (sample) {
      # same as Bestie:::SampleParameters
      beta_by_cluster[[k]] <- as.vector(
        mvtnorm::rmvt(
          1,
          sigma = sigma,
          df = df,
          delta = as.vector(mb)
        )
      )
    } else {
      beta_by_cluster[[k]] <- as.vector(mb)
    }
  }

  if (sample) {
    k <- sample(seq_len(K), 1, prob = scoreparam$clusterWeights)
    beta_by_cluster[[k]]
  } else {
    Reduce(
      "+",
      Map(
        function(beta, w) w * beta,
        beta_by_cluster,
        scoreparam$clusterWeights
      )
    )
  }
}

DP_DAGintervention <- function(incidences, dataParams, sample = TRUE) {
  one_dag <- function(incidence) {
    p <- ncol(incidence)
    coeffMatrix <- matrix(0, p, p)
    colnames(coeffMatrix) <- rownames(coeffMatrix) <- colnames(incidence)
    for (j in seq_len(p)) {
      parentNodes <- which(incidence[, j] == 1)
      if (length(parentNodes) == 0) next
      beta <- DP_sample_node_beta(
        j = j,
        parentNodes = parentNodes,
        dataParams = dataParams,
        sample = sample
      )
      coeffMatrix[parentNodes, j] <- beta
    }
    solve(diag(p) - coeffMatrix)
  }
  if (is.matrix(incidences)) {
    return(one_dag(incidences))
  }
  lapply(incidences, one_dag)
}
