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
  for (l in  1:L){
    clusterParams = clusterParams_sample[[l]]
    pointsPerCluster <- pointsPerCluster_sample[[l]]
    numLabels <- length(pointsPerCluster)
    probs <- matrix(0, nrow = N, ncol = numLabels)
    for (i in seq_len(N)) {
      probs[i, 1:numLabels] <- pointsPerCluster * 
        Likelihood(mdObj, 
                   y[i,, drop = FALSE], 
                   clusterParams)
    }
    probs_list[[l]] <- probs / rowSums(probs)
  }
  return(probs_list)
}
add_membershipp <- function(membershipp_list, membershipp, child, parents, active=T) {
  membershipp_list[[length(membershipp_list) + 1]] <- list(
    membershipp = membershipp,
    child = child, 
    parents = parents, 
    vars = c(child, parents),
    active = active
  )
  membershipp_list
}
adj_union_from_aliases <- function(aliases,  var_names) {
  n = length(var_names)
  adj <- matrix(0L, n, n, dimnames = list(var_names, var_names))
  for (i in seq_len(n)) {
    ai <- aliases[[i]]
    if (is.null(ai)) next
    if (!is.matrix(ai)) ai <- as.matrix(ai)
    
    # alias ids referenced for this node (drop NA/0 padding)
    parents <- unique(as.integer(ai[!is.na(ai) & ai != 0]))
    if (!length(parents)) next
    
    # Safety check: are these ids plausible row indices
    if (max(parents) > n) {
      stop(sprintf(
        "For node %d: max alias id (%d) > n (%d).",
        i, max(parents), i, n
      ))
    }
    
    
    if (length(parents)) adj[parents, i] <- 1L
  }
  
  adj
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
        vars    = dp_membershipp_list[[d]]$vars,
        active  = dp_membershipp_list[[d]]$active
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
  dp_scoreparam_list = Filter(function(e) all((param$labels[needed] %in% e$meta$vars) & e$meta$active), 
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
update_score_param <- function(param, space){
  for (i in 1:length(param$dp_scoreparam_list)){
    # deactivate all
    param$dp_scoreparam_list[[i]]$meta$active = F
  }
  # initiate new probability matrices
  Gamma_list <- list()
  for (child in colnames(space)){
    parents <- names(which(space[ , child] == 1))
    idx <- which(vapply(param$dp_scoreparam_list, function(e) {
      identical(e$meta$child, child) &&
        setequal(e$meta$parents, parents)
    }, logical(1)))
    
    # activate
    if (length(idx) > 0) {
      if (length(idx) > 1) stop("multiple DP with same child parents")
      for (j in idx) {
        param$dp_scoreparam_list[[j]]$meta$active <- TRUE
      }
    }
    
    # create new DP
    if (length(idx) == 0) {
      dp_data = param$data[,c(child, parents)]
      if (length(parents) == 0){
        dp <-  DirichletProcessGaussian(dp_data,
                                        alphaPriors = c(1, 20))
      }
      else{
        n_col = ncol(dp_data)
        g0Priors <- list(
          mu0    = rep(0, n_col),
          Lambda = diag(n_col) / param$T0scale,   # T = (1/t) I
          kappa0 = param$am,
          nu     = param$aw
        )
        
        dp <-  DirichletProcessMvnormal(dp_data, g0Priors)
      }
      dp <- Fit(dp, param$dp_iter)
      
      Gamma_sample <- dp_membership_probs(dp, param$dp_burnin, param$dp_n_sample)
      # add meta data in list
      Gamma_list <- add_membershipp(Gamma_list, 
                                    Gamma_sample, 
                                    child=child, 
                                    parents=parents, 
                                    active=TRUE)
      
    }
    
  }
  new_score = BiDAG::scoreparameters(scoretype = "usr", 
                                     data = param$data, 
                                     usrpar = list(pctesttype = "bge",
                                                   membershipp_list = Gamma_list,
                                                   am = param$am, 
                                                   aw = param$aw, 
                                                   T0scale = param$T0scale,
                                                   edgepf = 1
                                     )
  )
  # add new computed params to existing param
  param$dp_scoreparam_list <- append(param$dp_scoreparam_list, new_score$dp_scoreparam_list)
  param
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
## Sample DAGs with BiDAG
## ---------------------------------------------------------------------------------------------------------------------------------------------------------
sampleDAGs <- function(inData, nDigraphs = 50, seed=101, dname="", ...){
  n <- ncol(inData) # number of variables
  stepsave <- 10*round(n*n*log(n)) ## thinning interval, 
  ## to improve the properties of the chain (sample independence)
  iterations <- nDigraphs*stepsave ## nDAGs is the desired size of the DAGs ensemble
  
  if(!dir.exists("./saveout")) {dir.create("./saveout")} ## create output folder if it doesn't exist
  if(!file.exists(paste0("./saveout/dagdraw", n, "seed", seed, dname, ".RData"))){
    # initialise the score object needed for later functions
    scoreObject <- BiDAG::scoreparameters(data=inData, ..., 
                                   nodeslabels = colnames(inData))
    set.seed(seed) 
    ## set the seed for the generation of random numbers (for reproducibility)
    
    # find the search space with iterative search
    itFit <- BiDAG::iterativeMCMC(scoreObject, 
                           startorder = sample(c(1:n)), 
                           scoreout = TRUE) ## find iterative search space
    searchSpace <- itFit$endspace
    
    # sample a starting DAG with order MCMC
    # the default length of the chain is 6*n^2/log(n)
    # orderSample$score is the score for the highest scoring DAG found
    orderSample <- BiDAG::orderMCMC(scoreObject, 
                             MAP = FALSE, 
                             startspace = searchSpace, 
                             startorder = sample(c(1:n)), 
                             chainout = TRUE)
    startDAG <- dplyr::last(orderSample$traceadd$incidence) # extract selected (last) DAG
    startDAGscore <- dplyr::last(orderSample$trace) # extract score for the selected DAG
    
    # sample an ensemble of nDAGs DAGs from the posterior by using partition MCMC
    partitionSample <- BiDAG::partitionMCMC(scoreObject, 
                                     startspace = searchSpace, 
                                     startDAG = startDAG, 
                                     iterations = iterations, 
                                     stepsave = stepsave)
    
    ## extract the collection of sampled DAGs
    sampledDAGs <- partitionSample$traceadd$incidence 
    DAGscores <-partitionSample$trace # extract their scores
    
    ## save the the collection of sampled DAGs and their scores to an .RData file
    save(sampledDAGs, DAGscores, scoreObject,
         file=paste0("./saveout/dagdraw", n, "seed", seed, dname, ".RData"))
  }
}  

## ---------------------------------------------------------------------------------------------------------------------------------------------------------
## Internal function to set-up plot style
## Define some default plotting settings for the effects
## ---------------------------------------------------------------------------------------------------------------------------------------------------------
setPlot <- function(dt = NULL, ...){
  plot(dt, 
       axes=FALSE, 
       frame.plot=FALSE, 
       xlab="", 
       ylab="", 
       main="",
       ...)
}


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
## Plot a summary DAG from a collection of samples
## ---------------------------------------------------------------------------------------------------------------------------------------------------------
dagviz <- function(dags4plot, # roundprobs,
                   daglabs = colnames(dags4plot[[1]]), # character vector of names
                   nn = length(daglabs),
                   grouped_nodes = NULL, # list(sample(order(daglabs), 2)),
                   # each element in the above list 
                   # will have its nodes placed inside a group
                   ## It won't work with fewer than 2 nodes
                   node_colour = "#1e90ff", # node colour
                   font_colour = "#fffaf0", # font colour in nodes
                   edge_colour = "#68228b", # edge colour
                   edge_threshold = 0.1, # only show edges above this threshold
                   group_colour = "#add8e6", # colour for group of nodes
                   edge_width = 1.5, # edge width
                   title_text = "Directed Acyclic Graph\n",
                   # title text (use \n for line breaks)
                   style_mat = matrix(c(1,2), nrow=2*nn + 1, ncol=2*nn)[1:nn, 1:nn],
                   # checker-board for illustration
                   edge_styles = c("solid", "dashed")
                   ){
  ## calculate plotting parameters from the DAG data
  dagsarray <- round(simplify2array(lapply(dags4plot, as.matrix)), 8)
  edgeprobs <- apply(dagsarray, c(1,2), mean)
  roundprobs <- round(255*edgeprobs, 0) # rounded to 0-255 for line intensity
  
  ## Build graphviz code
  graphcommand <- paste0("digraph G { \n node [color=\"", 
                         node_colour, 
                         "\", style=filled, fontcolor=\"", 
                         font_colour, "\"]; \n ")
  if (length(grouped_nodes) > 0) {
    for (kk in 1:length(grouped_nodes)) {
      node_text1 <- paste0(daglabs[grouped_nodes[[kk]]], 
                           ";", collapse = " ")
      node_text2 <- paste0(daglabs[grouped_nodes[[kk]]], 
                           collapse = " ")
      graphcommand <- paste0(graphcommand,
                           "subgraph cluster", 
                           kk, 
                           " { \n bgcolor=\"", 
                           group_colour, 
                           "\" \n penwidth=0 \n", 
                           node_text1,
                           " \n {rank=same ", 
                           node_text2,"} \n } \n ") 
      }
    }
  for(ii in 1:nn){
    for(jj in 1:nn){
    if(edgeprobs[ii,jj] > edge_threshold){
      graphcommand<-paste0(graphcommand, 
                          daglabs[ii],
                          " -> ", 
                          daglabs[jj], 
                          " [color=\"", 
                          edge_colour, 
                          as.hexmode(roundprobs[ii,jj]),
                          "\", penwidth=", 
                          edge_width, 
                          ", style=", 
                          edge_styles[style_mat[ii, jj]],"] \n")
      }
    }
  }
  grViz(paste0(graphcommand,"labelloc=\"t\";\n label= \"", title_text, "\";}"))
}

## filename for gv file if desired, not needed with DiagrammeR
## gv_output <- "DAG.gv"
## write(graphcommand, file=gv_output)

## ---------------------------------------------------------------------------------------------------------------------------------------------------------
## grViz() by itself produces a html output not visible in pdfs (even by setting always_allow_html: true in the YAML)
## ---------------------------------------------------------------------------------------------------------------------------------------------------------
## Save summary DAG to a file and display on pdf
## ---------------------------------------------------------------------------------------------------------------------------------------------------------

displayDAG <- function(g2plot, figname = "dataDAG.png"){
  g2plot %>%
    export_svg %>%
    charToRaw %>%
    rsvg_png(figname, width = 12000, height = 12000)
  ## Increasing the size also increases the quality of the graph
  knitr::include_graphics(figname)
}
