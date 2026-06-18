## ---------------------------------------------------------------------------------------------------------------------------------------------------------
## Sample DAGs with BiDAG
## ---------------------------------------------------------------------------------------------------------------------------------------------------------
sampleDAGs <- function(inData, searchspace, weighted = FALSE, nDigraphs = 50, seed=101, dname="", ...){
  scoreObject <- searchspace$score
  scoreObject$data = inData
  scoreObject$nodeslabels = colnames(inData)

  n <- ncol(inData) # number of variables
  stepsave <- 10*round(n*n*log(n)) ## thinning interval, 
  ## to improve the properties of the chain (sample independence)
  iterations <- nDigraphs*stepsave ## nDAGs is the desired size of the DAGs ensemble
  
  if(!dir.exists("./Sachs/saveout")) {dir.create("./Sachs/saveout")} ## create output folder if it doesn't exist
  if(!file.exists(paste0("./Sachs/saveout/dagdraw", n, "seed", seed, dname, ".RData"))){
    # initialise the score object needed for later functions

    set.seed(seed) 
    ## set the seed for the generation of random numbers (for reproducibility)
    
    # find the search space with iterative search
    endspace <- searchspace$endspace
    
    # sample a starting DAG with order MCMC
    # the default length of the chain is 6*n^2/log(n)
    # orderSample$score is the score for the highest scoring DAG found
    orderSample <- orderMCMC(scoreObject, 
                             MAP = FALSE, 
                             startspace = endspace, 
                             chainout = TRUE,
                             compress = FALSE)
    startDAG <- last(orderSample$traceadd$incidence) # extract selected (last) DAG
    startDAGscore <- last(orderSample$trace) # extract score for the selected DAG
    
    # sample an ensemble of nDAGs DAGs from the posterior by using partition MCMC
    partitionSample <- partitionMCMC(scoreObject, 
                                     startspace = endspace, 
                                     startDAG = startDAG, 
                                     iterations = iterations, 
                                     stepsave = stepsave,
                                     compress = FALSE)
    
    ## extract the collection of sampled DAGs
    sampledDAGs <- partitionSample$traceadd$incidence 
    DAGscores <- partitionSample$trace # extract their scores
    if(weighted) {
      ndags <- length(partitionSample$trace)
      score_cache <- new.env(parent = emptyenv())
      dag_key <- function(dag) {
        paste(as.integer(dag), collapse = "")
      }
      weights <- numeric(ndags)
      target_scores <- numeric(ndags)
      
      for (k in seq_len(ndags)) {
        dag <- partitionSample$traceadd$incidence[[k]]
        key <- dag_key(dag)
        
        if (exists(key, envir = score_cache, inherits = FALSE)) {
          target_score <- get(key, envir = score_cache)
        } else {
          target_score <- DPscoreDAG(
            param = scoreObject,
            dag = dag
          )
          assign(key, target_score, envir = score_cache)
        }
        
        proposal_score <- partitionSample$trace[k]
        
        target_scores[k] <- target_score
        weights[k] <- target_score - proposal_score
      }
      partitionSample$weights <- weights - logSumExp(weights)  # normalize weights
      ind <- sample(seq_along(sampledDAGs), size = nDigraphs, replace = TRUE, prob = exp(partitionSample$weights))
      sampledDAGs = sampledDAGs[ind]
      DAGscores <- target_scores[ind]
    }
    
    ## save the the collection of sampled DAGs and their scores to an .RData file
    save(sampledDAGs, DAGscores, scoreObject,
         file=paste0("./Sachs/saveout/dagdraw", n, "seed", seed, dname, ".RData"))
  }
}  
makeTrueDAGSamples <- function(truegraph, nDAGs) {
  sampledDAGs <- replicate(nDAGs + 1, truegraph, simplify = FALSE)
  sampledDAGs <- lapply(sampledDAGs, function(A) {
    A <- as.matrix(A)
    colnames(A) <- rownames(A) <- colnames(truegraph)
    A
  })
  
  sampledDAGs
}
## ---------------------------------------------------------------------------------------------------------------------------------------------------------
## Compute intervention effects with Bestie
## ---------------------------------------------------------------------------------------------------------------------------------------------------------
computeEffects <- function(n, seed=101, dname="", DP=FALSE){
  
  if (!dir.exists("./Sachs/saveout")) {dir.create("./Sachs/saveout")} ## create output folder if it doesn't exist
  if(file.exists(paste0("./Sachs/saveout/dagdraw", n, "seed", seed, dname, ".RData"))){
    # load the collection of sampled DAGs
    load(file = paste0("./Sachs/saveout/dagdraw", n, "seed", seed, dname, ".RData"))

    # sample parameters and derive a matrix of intervention effects for each DAG, 
    # saving them all in a list;
    # First check if effects have already been estimated and saved to file
    if(!file.exists(paste0("./Sachs/saveout/effects", n, "seed", seed, dname, ".RData"))){
      if(DP){
        causalMats <- DP_DAGintervention(
          incidences = sampledDAGs,
          dataParams = scoreObject,
          sample = TRUE
        )
      } else{
        causalMats <- DAGintervention(sampledDAGs, scoreObject, sample=TRUE)
      }

      # save estimated intervention effects to an .RData file
      save(causalMats,
           file=paste0("./Sachs/saveout/effects", n, "seed", seed, dname, ".RData"))
    }
  } else(warning("Looking for causal effects without a DAG? Try runif()!"))
}
computeEffects_mem <- function(sampledDAGs, scoreObject, DP = FALSE) {
  if (DP) {
    DP_DAGintervention(
      incidences = sampledDAGs,
      dataParams = scoreObject,
      sample = TRUE
    )
  } else {
    DAGintervention(sampledDAGs, scoreObject, sample = TRUE)
  }
}
sampleTrueEffects <- function(Eff1, Eff2, n1, n2, labels = NULL, nSamples = 1000) {
  n <- nrow(Eff1)
  ntot <- n1 + n2
  out <- vector("list", nSamples)
  for (s in seq_len(nSamples)) {
    M <- if (runif(1) < n1 / ntot) Eff1 else Eff2
    
    if (!is.null(labels)) {
      colnames(M) <- rownames(M) <- labels
    }
    out[[s]] <- M
  }
  out
}
## ---------------------------------------------------------------------------------------------------------------------------------------------------------
## Load collection of sampled DAGs and effects
## ---------------------------------------------------------------------------------------------------------------------------------------------------------
loadsamples <- function(seeds, nn, dname = "") {
  nSeeds <- length(seeds)
  alldigraphs <- vector("list", nDAGs * nSeeds) # to store the graphs
  alleffs <- vector("list", nDAGs * nSeeds) # to store the matrices of effects
  for (nlevel in 1:nSeeds) {
    seednumber <- seeds[nlevel]
    ## Retrieve sampled DAGs - DAG chain
    load(file = paste0("./Sachs/saveout/dagdraw", nn, "seed", seednumber, dname, ".RData"))
    alldigraphs[1:nDAGs + (nlevel - 1) * nDAGs] <- sampledDAGs # remove the starting point
    ## Retrieve estimated effects
    load(file = paste0("./Sachs/saveout/effects", nn, "seed", seednumber, dname, ".RData"))
    alleffs[1:nDAGs + (nlevel - 1) * nDAGs] <- causalMats # remove the starting point
  }
  list(alldigraphs=alldigraphs, alleffs=alleffs)
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
                   edge_styles = c("solid", "dashed"),
                   rm_nodes = NULL
){
  ## calculate plotting parameters from the DAG data
  dagsarray <- round(simplify2array(dags4plot), 8)
  edgeprobs <- apply(dagsarray, c(1,2), mean)
  
  if (!is.null(rm_nodes)) {
    edgeprobs <- edgeprobs[-rm_nodes, -rm_nodes]
    nn <- ncol(edgeprobs)
    daglabs <- colnames(edgeprobs)
  }
  
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

## ---------------------------------------------------------------------------------------------------------------------------------------------------------
## Plot a matrix of posterior distributions of effects
## ---------------------------------------------------------------------------------------------------------------------------------------------------------
plotEffects <- function(effects4plot, 
                        sortlabs, ## an order vector of indiced
                        ## would it be better if names are saved in the list of effect matrices?
                        ## then use match if given a vector of order labels to find the indices in the Effect matrices
                        sigcutoff = 0.025,
                        xmargs = rep(0.1, 2),
                        label_size = 1,
                        title_text = "Distributions of Downstream Causal Effects\n"){
  orderingy = c(0, sortlabs) 
  ## Provide labels in sortlabs in a customised order if desired
  ### Pick an ordering which matches to some extent the node hierarchy on the DAG
  
  nn <- length(sortlabs)
  effsarray <- round(simplify2array(effects4plot), 8)
  aveffs <- apply(effsarray, c(1,2), mean)
  # estimate of the posterior mean of the intervention effects
  roundeffs <- round(aveffs, 3)
  
  efflabs <- colnames(roundeffs)
  
  minx <- min(apply(effsarray, c(1,2), min))
  maxx <- max(apply(effsarray, c(1,2), max))
  
  par(mfrow = c(nn+1, nn+1))
  par(oma = rep(0.2,4) + c(0,0,5*(title_text!=""),0))
  par(mar = rep(0,4))
  
  for(i in 0:nn+1){
    ii <- orderingy[i+1]
    for(j in 0:nn){
      jj <- orderingy[j+1]
      if(i==(nn+1) || j==0){
        setPlot(xlim=c(-1.1,1.1), ylim=c(0,1), col="darkorchid4")
        if(j==0 && i<(nn+1)){
          text(0,0.5, efflabs[ii], cex=label_size)
          }
        if(i==(nn+1) && j>0){
          text(0,0.5, efflabs[jj], srt=90, cex=label_size)
          }
        } else{
          if(!(roundeffs[ii,jj]==0 || roundeffs[ii,jj]==1) ){
            d <- density(effsarray[ii,jj,])
            setPlot(d, xlim=c(-2,2), col="dodgerblue")
            testquantiles <- quantile(effsarray[ii,jj,], c(sigcutoff, 1-sigcutoff))
            if(testquantiles[1]>0 || testquantiles[2]<0){
              u <- par("usr") # The coordinates of the plot area
              rect(u[1], u[3], u[2], u[4], col="#ffa07a44", border=NA)
              }
            polygon(d, col="dodgerblue", border="dodgerblue")
            abline(v=0, col="firebrick3", lty=2)
            u <- par("usr") # The coordinates of the plot area
            text( (u[1]+u[2])/2,(u[3]+u[4])/2, 
                  format(round(roundeffs[ii,jj], 2), nsmall = 2) )
            } else {
              setPlot(xlim=c(-2,2), ylim=c(0,1), col="darkorchid4")
              text(0, 0.5, roundeffs[ii,jj])
            }
        }
    }
  }
  if (title_text != "") {
    mtext(title_text, outer = TRUE, cex = 0.9)
  }
}
## ---------------------------------------------------------------------------------------------------------------------------------------------------------
## Plot compare distributions of effects
## ---------------------------------------------------------------------------------------------------------------------------------------------------------
wasserstein_effect_avg <- function(effects4plot,
                                   trueEffects,
                                   effectGraph,
                                   sortlabs = seq_len(nrow(effectGraph))) {
  effsarray <- round(simplify2array(effects4plot), 8)
  truearray <- round(simplify2array(trueEffects), 8)

  nn <- length(sortlabs)

  Wdist <- matrix(NA_real_, nn, nn)
  rownames(Wdist) <- rownames(effectGraph)[sortlabs]
  colnames(Wdist) <- colnames(effectGraph)[sortlabs]

  for (ii in sortlabs) {
    for (jj in sortlabs) {
      Wdist[ii, jj] <- wasserstein1d(
        truearray[ii, jj, ],
        effsarray[ii, jj, ]
      )
    }
  }

  keep <- effectGraph == 1

  list(
    Wdist = Wdist,
    avgW = mean(Wdist[keep], na.rm = TRUE),
    nEffects = sum(keep, na.rm = TRUE)
  )
}
plotCompareEffects <- function(effects4plot,
                                   trueEffects,
                                   sortlabs,
                                   sigcutoff = 0.025,
                                   xmargs = rep(0.1, 2),
                                   label_size = 1,
                                   title_text = "Distributions of Downstream Causal Effects\n") {
  
  orderingy = c(0, sortlabs) 
  ## Provide labels in sortlabs in a customised order if desired
  ### Pick an ordering which matches to some extent the node hierarchy on the DAG
  
  nn <- length(sortlabs)
  effsarray <- round(simplify2array(effects4plot), 8)
  truearray <- round(simplify2array(trueEffects), 8)
  aveffs <- apply(effsarray, c(1,2), mean)
  # estimate of the posterior mean of the intervention effects
  roundeffs <- round(aveffs, 3)
  
  efflabs <- colnames(roundeffs)
  
  minx <- min(apply(effsarray, c(1,2), min))
  maxx <- max(apply(effsarray, c(1,2), max))
  
  Wdist <- matrix(NA_real_, nn, nn)
    
  for (ii in sortlabs) {
    for (jj in sortlabs) {
      Wdist[ii, jj] <- wasserstein1d(
        truearray[ii, jj, ],
        effsarray[ii, jj, ]
      )
    }
  }
  
  Wscaled <- pmin(Wdist, 1)
  heatcols <- colorRampPalette(c("white", "gold", "orange", "red"))(100)
  
  par(mfrow = c(nn+1, nn+1))
  par(oma = rep(0.2,4) + c(0,0,5*(title_text!=""),0))
  par(mar = rep(0,4))
  
  for(i in 0:nn+1){
    ii <- orderingy[i+1]
    for(j in 0:nn){
      jj <- orderingy[j+1]
      if (i == (nn + 1) || j == 0) {
        setPlot(xlim = c(-1.1, 1.1), ylim = c(0, 1), col = "darkorchid4")
        if (j == 0 && i < (nn + 1)) {
          text(0, 0.5, efflabs[ii], cex = label_size)
        }
        if (i == (nn + 1) && j > 0) {
          text(0, 0.5, efflabs[jj], srt = 90, cex = label_size)
        }
      } else {
        heat_id <- max(1, ceiling(Wscaled[ii, jj] * 100))
        heat_col <- adjustcolor(heatcols[heat_id], alpha.f = 0.45)

        vals <- effsarray[ii, jj, ]
        if (sd(vals) > 1e-8) {
          d <- density(vals)
          setPlot(d, xlim = c(-2, 2), col = "dodgerblue")
          u <- par("usr")
          rect(u[1], u[3], u[2], u[4],
              col = heat_col,
              border = NA)

          polygon(d, col = "dodgerblue", border = "dodgerblue")
          abline(v = 0, col = "firebrick3", lty = 2)
          text(
            (u[1] + u[2]) / 2,
            (u[3] + u[4]) / 2,
            paste0(
              "W=",
              format(round(Wdist[ii, jj], 2), nsmall = 2)
            )
          )

        } else {
          setPlot(
            xlim = c(-2, 2),
            ylim = c(0, 1),
            col = "darkorchid4"
          )
          u <- par("usr")
          rect(
            u[1], u[3], u[2], u[4],
            col = heat_col,
            border = NA
          )
          text(
            0,
            0.5,
            paste0(
              "W=",
              format(round(Wdist[ii, jj], 2), nsmall = 2)
            )
          )
        }
      }
    }
  }
  if (title_text != "") {
    mtext(title_text, outer = TRUE, cex = 0.9)
  }
  
  invisible(Wdist)
}