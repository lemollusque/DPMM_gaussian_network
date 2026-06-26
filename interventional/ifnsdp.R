#---------------------- membership prob matrix ----------------------
dp_membership_probs <- function(fit, Y, N, L = 500) {
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
#----------------------  BiDAG ----------------------------------
usrscoreparameters <- function(initparam, usrpar = list(Imat = NULL, pctesttype = "bge", am = 1, edgepmat = NULL, bgremove = TRUE)){
  n <- initparam$n
  Imat <- usrpar$Imat
  bgn <- ncol(Imat)
  colnames(Imat) <- paste0("i", 1:bgn)
  data <- initparam$data
  weightvector <- usrpar$weightvector
  
  exps <- mgcv::uniquecombs(Imat) # experimental conditions
  expsrows <- attr(exps, "index") # rows with each condition

  # dirichlet params
  prior <- usrpar$dp_prior
  mcmc <- usrpar$dp_mcmc
  L <- usrpar$dp_n_sample
  dp_fits <- usrpar$dp_fits

  # dp on full data
  initparamlocal <- initparam
  initparamlocal$data <- data
  initparamlocal$n <- ncol(data)  
  
  # prepare dirichlet gamma list
  Gamma_list <- list()
  for (f in seq_len(dp_fits)) {
    output <- list(out_param = TRUE, out_type = "FULL")  
    fit <- PYdensity(y = data, mcmc = mcmc, prior = prior, output = output)
    Gamma_sample <- dp_membership_probs(fit, data , nrow(data), L)
    Gamma_list[[f]] <- list(membershipp = Gamma_sample)
  }
  usrpar$membershipp_list = Gamma_list
  dp_score <- dpscoreparameters(initparamlocal, usrpar = usrpar)
 
  if (bgn > 0 && nrow(exps) > 1) {
    initparam <- scoreparameters(scoretype = "bge", data = cbind(Imat,data), bgnodes = 1:ncol(Imat), bgepar = list(am = usrpar$am),
                                 edgepmat = usrpar$edgepmat)
    initparam$dp_score <- dp_score
    initparam$type <- "usr" # make sure it knows that we have redefined the score
    initparam$pctesttype <- "bge"
    initparam$bgremove <- usrpar$bgremove

    make_key <- function(idx) paste(sort(idx), collapse = "_")
    all_subset_keys <- character(0)
    all_subsets <- list()
    for (r in seq_len(bgn)) {
      for (ip in combn(seq_len(bgn), r, simplify = FALSE)) {
        local_exps <- attr(mgcv::uniquecombs(exps[, ip, drop = FALSE]), "index")
        for (ii in seq_len(max(local_exps))) {
          idx <- which(local_exps == ii)
          key <- make_key(idx)
          if (!(key %in% all_subset_keys)) {
            all_subset_keys <- c(all_subset_keys, key)
            all_subsets[[key]] <- idx
          }
        }
      }
    }

    dp_scores <- list()
    for (key in names(all_subsets)) {
      idx_conditions <- all_subsets[[key]]
      rows <- which(expsrows %in% idx_conditions) 
      datalocal <- data[rows, , drop = FALSE]
      initparamlocal <- initparam
      initparamlocal$data <- datalocal
      initparamlocal$n <- ncol(datalocal)
      initparamlocal$N <- nrow(datalocal)
      # prepare dirichlet gamma list
      Gamma_list <- list()
      for (f in seq_len(dp_fits)) {
        output <- list(out_param = TRUE, out_type = "FULL")  
        fit <- PYdensity(y = datalocal, mcmc = mcmc, prior = prior, output = output)
        Gamma_sample <- dp_membership_probs(fit, datalocal , nrow(datalocal), L)
        Gamma_list[[f]] <- list(membershipp = Gamma_sample)
      }
      usrparlocal <- usrpar
      usrparlocal$membershipp_list <- Gamma_list
      
      dp_scores[[key]] <- list(
        idx_conditions = idx_conditions,
        score = dpscoreparameters(initparamlocal, usrparlocal)
      )
    }
    initparam$dp_scores <- dp_scores
  } else {
    initparam$dp_scores <- list(dp_score)
  }
  
  initparam$obsdata <- data
  initparam$usrpar <- usrpar
  initparam$exps <- exps
  initparam$expsrows <- expsrows
  initparam
}
dpscoreparameters <- function(initparam, 
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
  bgeinitparam <- BGeaugment(covmat, means, N, n, usrpar$am, usrpar$aw, initparam$logedgepmat)
  initparam$bgeinitparam <- bgeinitparam

  initparam
}
### This function evaluates the log score of a node given its parents
usrDAGcorescore <- function (j, parentnodes, n, param) {
  bgn = param$bgn
  iparents <- parentnodes[which(parentnodes %in% param$bgnodes)]
  if (length(iparents) == 0 || nrow(param$exps) < 2) {# use standard BGe score
    localparam <- param$dp_score$bgeinitparam
    outscore <- BiDAG:::DAGcorescore(j - bgn, parentnodes - bgn, n, localparam)
  } else {
    parents <- setdiff(parentnodes, iparents)
    # find the different exp conditions for these parents
    local_exps <- attr(mgcv::uniquecombs(param$exps[, iparents, drop = FALSE]), "index")
    outscore <- 0
    for (ii in 1:max(local_exps)) {
      key <- paste(sort(which(local_exps == ii)), collapse = "_")
      localparam <- param$dp_scores[[key]]$score$bgeinitparam
      if (length(parents) > 0) {
        outscore <- outscore + BiDAG:::DAGcorescore(1, 1:length(parents) + 1, param$n, localparam)
      }  else {
        outscore <- outscore + BiDAG:::DAGcorescore(1, parents, param$n, localparam)
      }
    }
  }
  outscore
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
  usrpar = param$usrpar
  bgn = param$bgn
  n <- ncol(dag)
  L <- usrpar$dp_n_sample
  dag_scores <- numeric(L)

  for (j in (bgn + 1):n) {
    parentnodes <- which(dag[, j] == 1)
    iparents <- parentnodes[which(parentnodes %in% param$bgnodes)]
    if (length(iparents) == 0 || nrow(param$exps) < 2) {# use standard BGe score
      localparam <- param$dp_score
      for (l in 1:L) {
        dag_scores[l] <- dag_scores[l] + targetCoreScore(
          j - bgn,
          parentnodes - bgn,
          n - bgn,
          localparam,
          l
        )
      }
    } else {
      parents <- setdiff(parentnodes, iparents)
      # find the different exp conditions for these parents
      local_exps <- attr(mgcv::uniquecombs(param$exps[, iparents, drop = FALSE]), "index")
      for (ii in 1:max(local_exps)) {
        key <- paste(sort(which(local_exps == ii)), collapse = "_")
        dpscore <- param$dp_scores[[key]]$score
        for (l in 1:L) {
          dag_scores[l] <- dag_scores[l] + targetCoreScore(
            j - bgn,
            parents - bgn,
            n - bgn,
            dpscore,
            l
          )
        }
      }
    }         
  }
  logMeanExpLogs(dag_scores)
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
# effect estimation ##############################################################
# Bestie
DP_sample_node_beta <- function(j, parentNodes, dataParams, sample = TRUE) {
  # Chose one DP draw for the scoring:
  # If several DP fits exist for the same child, choose one
  d <- sample(1:length(dataParams$dp_scoreparam_list), 1)
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
    usrpar <- dataParams$usrpar
    L <- usrpar$dp_n_sample
    bgn <- dataParams$bgn
    n <- ncol(incidence)
    p <- ncol(incidence)
    coeffMatrix <- matrix(0, p - bgn, p - bgn)
    colnames(coeffMatrix) <- rownames(coeffMatrix) <- colnames(incidence)[(bgn + 1):p]
    for (j in (bgn + 1):p) {
      parentNodes <- which(incidence[, j] == 1)
      iparents <- parentNodes[which(parentNodes %in% dataParams$bgnodes)]
      parents <- setdiff(parentNodes, iparents)
      if (length(parents) == 0) next
      if (length(iparents) == 0 || nrow(dataParams$exps) < 2) {# use standard BGe score
        localparam <- dataParams$dp_score
      } 
      else {
        # find the different exp conditions for these parents
        exp_conds <- mgcv::uniquecombs(dataParams$exps[, iparents, drop = FALSE])
        local_exps <- attr(exp_conds, "index")
        ii <- which(rowSums(exp_conds) == 0) # this defines the observational state
        key <- paste(sort(which(local_exps == ii)), collapse = "_")
        localparam <- dataParams$dp_scores[[key]]$score
      }

      beta <- DP_sample_node_beta(
        j = j - bgn,
        parentNodes = parents - bgn,
        dataParams = localparam,
        sample = sample
      )
      coeffMatrix[parents - bgn, j - bgn] <- beta
    }
    solve(diag(p - bgn) - coeffMatrix)
  }
  if (is.matrix(incidences)) {
    return(one_dag(incidences))
  }
  lapply(incidences, one_dag)
}
