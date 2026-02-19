#dirichletprocess
# dirichletprocess:::Fit.default
dirichletprocess_Fit_default = function (dpObj, its, updatePrior = FALSE, progressBar = interactive()) 
{
  if (progressBar) {
    pb <- txtProgressBar(min = 0, max = its, width = 50, 
                         char = "-", style = 3)
  }
  alphaChain <- numeric(its)
  likelihoodChain <- numeric(its)
  weightsChain <- vector("list", length = its)
  clusterParametersChain <- vector("list", length = its)
  priorParametersChain <- vector("list", length = its)
  labelsChain <- vector("list", length = its)
  for (i in seq_len(its)) {
    alphaChain[i] <- dpObj$alpha
    weightsChain[[i]] <- dpObj$pointsPerCluster/dpObj$n
    clusterParametersChain[[i]] <- dpObj$clusterParameters
    priorParametersChain[[i]] <- dpObj$mixingDistribution$priorParameters
    labelsChain[[i]] <- dpObj$clusterLabels
    likelihoodChain[i] <- sum(log(LikelihoodDP(dpObj)))
    dpObj <- ClusterComponentUpdate(dpObj)
    dpObj <- ClusterParameterUpdate(dpObj)
    dpObj <- UpdateAlpha(dpObj)
    if (updatePrior) {
      dpObj$mixingDistribution <- PriorParametersUpdate(dpObj$mixingDistribution, 
                                                        dpObj$clusterParameters)
    }
    if (progressBar) {
      setTxtProgressBar(pb, i)
    }
  }
  dpObj$weights <- dpObj$pointsPerCluster/dpObj$n
  dpObj$alphaChain <- alphaChain
  dpObj$likelihoodChain <- likelihoodChain
  dpObj$weightsChain <- weightsChain
  dpObj$clusterParametersChain <- clusterParametersChain
  dpObj$priorParametersChain <- priorParametersChain
  dpObj$labelsChain <- labelsChain
  if (progressBar) {
    close(pb)
  }
  return(dpObj)
}

# getS3method("ClusterComponentUpdate", "conjugate")
ClusterComponentUpdate_conjugate <-function (dpObj) 
{
    y <- dpObj$data
    n <- dpObj$n
    alpha <- dpObj$alpha
    clusterLabels <- dpObj$clusterLabels
    clusterParams <- dpObj$clusterParameters
    numLabels <- dpObj$numberClusters
    mdObj <- dpObj$mixingDistribution
    pointsPerCluster <- dpObj$pointsPerCluster
    predictiveArray <- dpObj$predictiveArray
    for (i in seq_len(n)) {
        currentLabel <- clusterLabels[i]
        pointsPerCluster[currentLabel] <- pointsPerCluster[currentLabel] - 
            1
        probs <- c(pointsPerCluster * Likelihood(mdObj, y[i, 
            , drop = FALSE], clusterParams), alpha * predictiveArray[i])
        probs[is.na(probs)] <- 0
        if (all(probs == 0)) {
            probs <- rep_len(1, length(probs))
        }
        newLabel <- sample.int(numLabels + 1, 1, prob = probs)
        dpObj$pointsPerCluster <- pointsPerCluster
        dpObj <- ClusterLabelChange(dpObj, i, newLabel, currentLabel)
        pointsPerCluster <- dpObj$pointsPerCluster
        clusterLabels <- dpObj$clusterLabels
        clusterParams <- dpObj$clusterParameters
        numLabels <- dpObj$numberClusters
    }
    dpObj$pointsPerCluster <- pointsPerCluster
    dpObj$clusterLabels <- clusterLabels
    dpObj$clusterParameters <- clusterParams
    dpObj$numberClusters <- numLabels
    return(dpObj)
}