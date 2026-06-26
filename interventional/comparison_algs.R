# Compare the fit to the true graph
compareFit <- function(graphlist, truegraph, weights = NULL, bgn=c(0)) {
  graphlist <- lapply(graphlist, function(x) x[-bgn, -bgn])
  cpdagList <- NULL  # check if DAGs are acyclic
  try(cpdagList <- lapply(graphlist, function(x) BiDAG:::dagadj2cpadj(x)), silent = T)
  pattList <- lapply(cpdagList, function(x) pdag2pattern(x))
  trueCPDAG <- BiDAG:::dagadj2cpadj(truegraph)
  truepatt <- pdag2pattern(trueCPDAG)

  SHD <- vector()
  TP <- vector()
  FP <- vector()
  
  for(pat in pattList) {  # Compute expected SHD, TP & FP
    comp <- compareGs(as.matrix(pat), truepatt)
    SHD <- append(SHD, comp["SHD"])
    TP <- append(TP, comp["TP"])
    FP <- append(FP, comp["FP"])
  }
  
  if(!is.null(weights)) {
    eshd <- sum(SHD * exp(weights))
    eTP <- sum(TP * exp(weights))
    eFP <- sum(FP * exp(weights))
  }
  
  else {
    eshd <- mean(SHD)
    eTP <- mean(TP)
    eFP <- mean(FP)
  }
  c(eshd, eTP, eFP)
}

# This function extracts the skeleton from a graph
Gskel <- function(incidence) {
  1*(incidence|t(incidence))
}

# This function compares an estimated graph to the true one
compareGs <- function (estG, trueG) {
  estSkel <- Gskel(estG) # estimated skeleton
  trueSkel <- Gskel(trueG) # true skeleton
  P <- sum(trueSkel)/2 # number of positives
  diffSkel <- estSkel - trueSkel
  extra_edges <- which(diffSkel > 0) # edges in estimated but not true EG
  FP <- length(extra_edges)/2 # count to FPs
  estG[extra_edges] <- 0 # remove them from further comparisons
  missing_edges <- which(diffSkel < 0) # edges in true but not estimated EG
  FN <- length(missing_edges)/2 # count to FNs
  trueG[missing_edges] <- 0 # remove them from further comparisons
  # modified graphs have the same skeletons, so now just need to count mismatches
  mismatches <- 1*(estG != trueG)
  wrong_order <- sum(Gskel(mismatches))/2 # number of wrongly oriented edges
  FP <- FP + wrong_order/2 # include half in FP
  FN <- FN + wrong_order/2 # and half in FN
  SHD <- FP + FN # shd is the sum of errors
  TP <- P - FN # true positives are without false negatives
  # TPR, FPR_P
  if (P == 0) { # true graph is empty
    if (FP >= 0) {
      TPR <- 0
      FPR_P <- 1
    } else {
      TPR <- 1
      FPR_P <- 0
    }
  } else { # true graph is non-empty
    TPR <- TP/P
    FPR_P <- FP/P
  }
  compGs <- c(TP, FP, SHD, TPR, FPR_P, P)
  names(compGs) <- c("TP","FP", "SHD", "TPR", "FPR_P", "P")
  return(compGs)
}
# Convert a cpdag to a pattern graph
pdag2pattern <- function(G) {
  pattern_graph <- matrix(0, nrow(G), ncol(G))
  ind <- which(G == 1, arr.ind = TRUE)
  
  if(length(ind) != 0) {
    for (i in 1:nrow(ind)) {
      x <- ind[i, 1]
      y <- ind[i, 2]
      allZ <- setdiff(which(G[ ,y] == 1), x)  # nodes directed towards y excluding x
      for (z in allZ) {
        if(G[y,x] == 0 && G[y,z] == 0 && G[x,z] == 0 && G[z,x] == 0) {
          pattern_graph[x,y] <- pattern_graph[z,y] <- 1  # save v-structure x -> y <- z
          G[x,y] <- G[z,y] <- 0  # delete v-structure from old pdag
        }
      }
    }
  }
  pdag <- 1 * (G | t(G))
  return(pattern_graph + pdag)
}
