#bidag fns
bidag_scoreparameters = function (scoretype = c("bge", "bde", "bdecat", "usr"), data, 
                      bgepar = list(am = 1, aw = NULL, edgepf = 1), bdepar = list(chi = 0.5, 
                                                                                  edgepf = 2), bdecatpar = list(chi = 0.5, edgepf = 2), 
                      dbnpar = list(samestruct = TRUE, slices = 2, b = 0, stationary = TRUE, 
                                    rowids = NULL, datalist = NULL, learninit = TRUE), usrpar = list(pctesttype = c("bge", 
                                                                                                                    "bde", "bdecat")), mixedpar = list(nbin = 0), MDAG = FALSE, 
                      DBN = FALSE, weightvector = NULL, bgnodes = NULL, edgepmat = NULL, 
                      nodeslabels = NULL) 
{
  initparam <- list()
  if (DBN) {
    dbnpardef <- list(samestruct = TRUE, slices = 2, b = 0, 
                      stationary = TRUE, rowids = NULL, datalist = NULL, 
                      learninit = TRUE)
    dbnpardef[names(dbnpar)] <- dbnpar[names(dbnpar)]
    dbnpar <- dbnpardef
    if (is.null(dbnpar$b)) 
      bgnodes <- NULL
    else if (dbnpar$b > 0) 
      bgnodes <- c(1:dbnpar$b)
    else bgnodes <- NULL
    initparam$learninit <- dbnpar$learninit
    if (!is.null(dbnpar$samestruct)) {
      initparam$split <- !dbnpar$samestruct
    }
    else {
      initparam$split <- FALSE
    }
    if (dbnpar$slices > 2 & !dbnpar$stationary) 
      MDAG <- TRUE
  }
  bgn <- length(bgnodes)
  if (DBN) {
    if (is.null(dbnpar$datalist)) {
      n <- (ncol(data) - bgn)/dbnpar$slices + bgn
    }
    else {
      n <- (ncol(data[[2]]) - bgn)/2 + bgn
    }
  }
  else n <- ncol(data)
  nsmall <- n - bgn
  if (!(scoretype %in% c("bge", "bde", "bdecat", "usr", "mixed"))) {
    stop("Scoretype should be bge (for continuous data), bde (for binary data) bdecat (for categorical data) or usr (for user defined)")
  }
  if (!DBN) {
    if (anyNA(data)) {
      stop("Dataset contains missing data")
    }
    if (ncol(data) != nsmall + bgn) {
      stop("n and the number of columns in the data do not match")
    }
  }
  else {
    if (dbnpar$stationary) {
      if (is.null(dbnpar$datalist)) {
        if (ncol(data) != nsmall * dbnpar$slices + bgn) {
          stop("n, bgn and the number of columns in the data do not match")
        }
      }
    }
  }
  if (!is.null(weightvector)) {
    if (length(weightvector) != nrow(data)) {
      stop("Length of the weightvector does not match the number of rows (observations) in data")
    }
  }
  if (scoretype == "bde") {
    if (!all(sapply(data, function(x) x %in% c(0, 1, NA)))) {
      stop("Dataset contains non-binary values")
    }
  }
  if (scoretype == "bdecat") {
    indx <- sapply(data, is.factor)
    data[indx] <- lapply(data[indx], function(x) as.numeric(x) - 
                           1)
    if (!all(unlist(lapply(data, function(x) setequal(unique(x), 
                                                      c(0:max(x))))))) {
      stop("Some variable levels are not present in the data")
    }
  }
  if (is.null(nodeslabels)) {
    if (!DBN) {
      if (all(is.character(colnames(data)))) {
        nodeslabels <- colnames(data)
      }
      else {
        nodeslabels <- sapply(c(1:n), function(x) paste("v", 
                                                        x, sep = ""))
      }
    }
    else {
      if (dbnpar$stationary & is.null(dbnpar$datalist)) {
        if (all(is.character(colnames(data)))) {
          nodeslabels <- colnames(data)
        }
        else {
          if (!is.null(bgnodes)) {
            staticnames <- sapply(c(1:bgn), function(x) paste("s", 
                                                              x, sep = ""))
            dynamicnames <- rep(sapply(c(1:nsmall), function(x) paste("v", 
                                                                      x, sep = "")), dbnpar$slices)
            for (i in 2:dbnpar$slices) {
              dynamicnames[1:nsmall + (i - 1) * nsmall] <- paste(dynamicnames[1:nsmall + 
                                                                                (i - 1) * nsmall], ".", i, sep = "")
            }
            nodeslabels <- c(staticnames, dynamicnames)
          }
          else {
            nodeslabels <- rep(sapply(c(1:n), function(x) paste("v", 
                                                                x, sep = "")), dbnpar$slices)
            for (i in 2:dbnpar$slices) {
              nodeslabels[1:nsmall + (i - 1) * nsmall] <- paste(nodeslabels[1:nsmall + 
                                                                              (i - 1) * nsmall], ".", i, sep = "")
            }
          }
        }
      }
      else {
        nodeslabels <- colnames(data[[2]])
      }
    }
  }
  multwv <- NULL
  if (is.null(dbnpar$datalist)) 
    colnames(data) <- nodeslabels
  initparam$labels <- nodeslabels
  initparam$type <- scoretype
  initparam$DBN <- DBN
  initparam$MDAG <- MDAG
  initparam$weightvector <- weightvector
  initparam$data <- data
  if (DBN) {
    initparam$bgnodes <- c(1:n + nsmall)
    if (bgn > 0) {
      initparam$static <- c(1:bgn)
    }
    initparam$mainnodes <- c(1:nsmall + bgn)
  }
  else {
    initparam$bgnodes <- bgnodes
    initparam$static <- bgnodes
    if (!is.null(bgnodes)) {
      initparam$mainnodes <- c(1:n)[-bgnodes]
    }
    else initparam$mainnodes <- c(1:n)
  }
  initparam$bgn <- bgn
  initparam$n <- n
  initparam$nsmall <- nsmall
  if (DBN) {
    if (dbnpar$stationary) {
      initparam$labels.short <- initparam$labels[1:(n + 
                                                      nsmall)]
    }
    else {
      nodeslabels <- colnames(data[[1]])
      initparam$labels <- nodeslabels
      initparam$labels.short <- colnames(data[[1]])
    }
  }
  else {
    initparam$labels.short <- initparam$labels
  }
  if (is.null(edgepmat)) {
    initparam$logedgepmat <- NULL
  }
  else {
    if (all(edgepmat > 0)) {
      initparam$logedgepmat <- log(edgepmat)
    }
    else stop("all entries of edgepmat matrix must be bigger than 0! 1 corresponds to no penalization")
  }
  if (DBN) {
    if (!dbnpar$stationary) {
      initparam$stationary <- FALSE
      initparam$slices <- length(data)
      initparam$intstr <- list()
      initparam$trans <- list()
      initparam$usrinitstr <- list()
      initparam$usrintstr <- list()
      initparam$usrtrans <- list()
      initparam$usrinitstr$rows <- c(1:n)
      initparam$usrinitstr$cols <- c(1:nsmall + bgn)
      if (bgn == 0) 
        initparam$usrintstr$rows <- c(1:nsmall + n)
      else initparam$usrintstr$rows <- c(1:bgn, 1:nsmall + 
                                           n)
      initparam$usrintstr$cols <- c(1:nsmall + n)
      initparam$usrtrans$rows <- c(1:nsmall + bgn)
      initparam$usrtrans$cols <- c(1:nsmall + n)
      if (bgn != 0) {
        initparam$intstr$rows <- c(1:bgn + nsmall, 1:nsmall)
      }
      else {
        initparam$intstr$rows <- c(1:nsmall)
      }
      initparam$intstr$cols <- c(1:nsmall)
      initparam$trans$rows <- c(1:nsmall + n)
      initparam$trans$cols <- c(1:nsmall)
      initparam$paramsets <- list()
      if (!is.null(edgepmat)) {
        edgepmatfirst <- edgepmat[1:n, 1:n]
        edgepmat <- DBNbacktransform(edgepmat, initparam, 
                                     nozero = TRUE)$trans
      }
      else {
        edgepmatfirst <- NULL
      }
      initparam$nsets <- length(data)
      datalocal <- data[[length(data)]]
      if (bgn > 0) 
        datalocal <- datalocal[, c(1:nsmall + bgn, 1:bgn)]
      initparam$paramsets[[length(data)]] <- scoreparameters(scoretype = scoretype, 
                                                             datalocal, weightvector = NULL, bgnodes = NULL, 
                                                             bgepar = bgepar, bdepar = bdepar, bdecatpar = bdecatpar, 
                                                             dbnpar = dbnpar, edgepmat = edgepmatfirst, DBN = FALSE)
      for (i in 1:(length(data) - 1)) {
        datalocal <- data[[i]]
        if (bgn > 0) 
          datalocal <- datalocal[, c(1:nsmall + nsmall + 
                                       bgn, 1:bgn, 1:nsmall + bgn)]
        else {
          datalocal <- datalocal[, c(1:nsmall + nsmall, 
                                     1:nsmall)]
        }
        initparam$paramsets[[i]] <- scoreparameters(scoretype = scoretype, 
                                                    datalocal, weightvector = NULL, bgnodes = initparam$bgnodes, 
                                                    bgepar = bgepar, bdepar = bdepar, bdecatpar = bdecatpar, 
                                                    dbnpar = dbnpar, edgepmat = edgepmat, DBN = FALSE)
      }
    }
    else {
      initparam$stationary <- TRUE
      initparam$slices <- dbnpar$slices
      if (!is.null(dbnpar$datalist)) {
        datalocal <- data[[2]]
        collabels <- colnames(datalocal)
        if (bgn > 0) 
          newbgnodes <- bgnodes + nsmall
        else newbgnodes <- bgnodes
      }
      else {
        datalocal <- data[, 1:(2 * nsmall + bgn)]
        collabels <- colnames(datalocal)
        if (bgn > 0) {
          bgdata <- data[, bgnodes]
          if (dbnpar$slices > 2) {
            for (jj in 1:(dbnpar$slices - 2)) {
              datatobind <- cbind(bgdata, data[, nsmall * 
                                                 jj + 1:(2 * nsmall) + bgn])
              colnames(datatobind) <- collabels
              datalocal <- rbind(datalocal, datatobind)
            }
          }
          newbgnodes <- bgnodes + nsmall
        }
        else {
          if (dbnpar$slices > 2) {
            for (jj in 1:(dbnpar$slices - 2)) {
              datatobind <- data[, n * jj + 1:(2 * n)]
              colnames(datatobind) <- collabels
              datalocal <- rbind(datalocal, datatobind)
            }
          }
          newbgnodes <- bgnodes
        }
      }
      if (bgn > 0) {
        datalocal <- datalocal[, c(1:nsmall + nsmall + 
                                     bgn, 1:bgn, 1:nsmall + bgn)]
      }
      else {
        datalocal <- datalocal[, c(1:n + n, 1:n)]
      }
      initparam$intstr <- list()
      initparam$trans <- list()
      initparam$usrinitstr <- list()
      initparam$usrintstr <- list()
      initparam$usrtrans <- list()
      initparam$usrinitstr$rows <- c(1:n)
      initparam$usrinitstr$cols <- c(1:nsmall + bgn)
      if (bgn == 0) 
        initparam$usrintstr$rows <- c(1:nsmall + n)
      else initparam$usrintstr$rows <- c(1:bgn, 1:nsmall + 
                                           n)
      initparam$usrintstr$cols <- c(1:nsmall + n)
      initparam$usrtrans$rows <- c(1:nsmall + bgn)
      initparam$usrtrans$cols <- c(1:nsmall + n)
      if (bgn != 0) {
        initparam$intstr$rows <- c(1:bgn + nsmall, 1:nsmall)
      }
      else {
        initparam$intstr$rows <- c(1:nsmall)
      }
      initparam$intstr$cols <- c(1:nsmall)
      initparam$trans$rows <- c(1:nsmall + n)
      initparam$trans$cols <- c(1:nsmall)
      if (!is.null(weightvector)) {
        weightvector.other <- rep(weightvector, dbnpar$slices - 
                                    1)
      }
      else {
        weightvector.other <- weightvector
      }
      lNA <- 0
      if (anyNA(datalocal)) {
        NArows <- which(apply(datalocal, 1, anyNA) == 
                          TRUE)
        lNA <- length(NArows)
        datalocal <- datalocal[-NArows, ]
        if (!is.null(weightvector)) {
          weightvector.other <- weightvector.other[-NArows]
        }
      }
      if (!is.null(edgepmat)) {
        edgepmatfirst <- edgepmat[1:n, 1:n]
        edgepmat <- DBNbacktransform(edgepmat, initparam, 
                                     nozero = TRUE)
        if (initparam$split) {
          edgepmat <- edgepmat$trans
        }
        else {
          initparam$logedgepmat <- log(edgepmat)
        }
      }
      else {
        edgepmatfirst <- NULL
      }
      initparam$otherslices <- scoreparameters(scoretype = scoretype, 
                                               datalocal, weightvector = weightvector.other, 
                                               bgnodes = initparam$bgnodes, bgepar = bgepar, 
                                               bdepar = bdepar, bdecatpar = bdecatpar, dbnpar = dbnpar, 
                                               edgepmat = edgepmat, DBN = FALSE)
      bdecatpar$edgepf <- 1
      bdepar$edgepf <- 1
      if (!is.null(dbnpar$datalist)) {
        datalocal <- data[[1]]
      }
      else {
        datalocal <- data[, 1:(nsmall + bgn)]
      }
      if (bgn == 0) {
        datalocal <- datalocal[, c(1:nsmall)]
      }
      else {
        datalocal <- datalocal[, c(1:nsmall + bgn, 1:bgn)]
      }
      if (!is.null(dbnpar$rowids)) {
        datalocal <- datalocal[which(dbnpar$rowids == 
                                       1), ]
      }
      if (anyNA(datalocal)) {
        NArows <- which(apply(datalocal, 1, anyNA) == 
                          TRUE)
        lNA <- lNA + length(NArows)
        datalocal <- datalocal[-NArows, ]
        if (!is.null(weightvector)) {
          weightvector <- weightvector[-NArows]
        }
      }
      if (lNA > 0) {
        cat(paste(lNA, "rows were removed due to missing data"), 
            "\n")
      }
      initparam$firstslice <- scoreparameters(scoretype = scoretype, 
                                              datalocal, weightvector = weightvector, bgnodes = newbgnodes, 
                                              bgepar = bgepar, bdepar = bdepar, bdecatpar = bdecatpar, 
                                              dbnpar = dbnpar, edgepmat = edgepmatfirst, DBN = FALSE)
    }
  }
  else if (scoretype == "bge") {
    if (is.null(bgepar$am)) {
      bgepar$am <- 1
    }
    if (is.null(bgepar$aw)) {
      bgepar$aw <- n + bgepar$am + 1
    }
    if (is.null(bgepar$edgepf)) {
      bgepar$edgepf <- 1
    }
    if (is.null(weightvector)) {
      N <- nrow(data)
      covmat <- cov(data) * (N - 1)
      means <- colMeans(data)
    }
    else {
      N <- sum(weightvector)
      forcov <- cov.wt(data, wt = weightvector, cor = TRUE, 
                       method = "ML")
      covmat <- forcov$cov * N
      means <- forcov$center
    }
    initparam$am <- bgepar$am
    initparam$aw <- bgepar$aw
    initparam$pf <- bgepar$edgepf
    initparam$N <- N
    initparam$means <- means
    mu0 <- numeric(n)
    T0scale <- bgepar$am * (bgepar$aw - n - 1)/(bgepar$am + 
                                                  1)
    T0 <- diag(T0scale, n, n)
    initparam$TN <- T0 + covmat + ((bgepar$am * N)/(bgepar$am + 
                                                      N)) * (mu0 - means) %*% t(mu0 - means)
    initparam$awpN <- bgepar$aw + N
    constscorefact <- -(N/2) * log(pi) + (1/2) * log(bgepar$am/(bgepar$am + 
                                                                  N))
    initparam$muN <- (N * means + bgepar$am * mu0)/(N + bgepar$am)
    initparam$SigmaN <- initparam$TN/(initparam$awpN - n - 
                                        1)
    initparam$scoreconstvec <- numeric(n)
    for (j in (1:n)) {
      awp <- bgepar$aw - n + j
      initparam$scoreconstvec[j] <- constscorefact - lgamma(awp/2) + 
        lgamma((awp + N)/2) + ((awp + j - 1)/2) * log(T0scale) - 
        j * log(initparam$pf)
    }
  }
  else if (scoretype == "bde") {
    if (is.null(bdepar$chi)) {
      bdepar$chi <- 0.5
    }
    if (is.null(bdepar$edgepf)) {
      bdepar$edgepf <- 2
    }
    if (is.null(weightvector)) {
      initparam$N <- nrow(data)
      initparam$d1 <- data
      initparam$d0 <- (1 - data)
    }
    else {
      initparam$N <- sum(weightvector)
      initparam$d1 <- data * weightvector
      initparam$d0 <- (1 - data) * weightvector
    }
    maxparents <- n - 1
    initparam$scoreconstvec <- rep(0, maxparents + 1)
    initparam$chi <- bdepar$chi
    initparam$pf <- bdepar$edgepf
    for (i in 0:maxparents) {
      noparams <- 2^i
      initparam$scoreconstvec[i + 1] <- noparams * lgamma(initparam$chi/noparams) - 
        2 * noparams * lgamma(initparam$chi/(2 * noparams)) - 
        i * log(initparam$pf)
    }
  }
  else if (scoretype == "bdecat") {
    if (is.null(bdecatpar$chi)) {
      bdecatpar$chi <- 0.5
    }
    if (is.null(bdecatpar$edgepf)) {
      bdecatpar$edgepf <- 2
    }
    maxparents <- n - 1
    initparam$chi <- bdecatpar$chi
    initparam$pf <- bdecatpar$edgepf
    initparam$scoreconstvec <- -c(0:maxparents) * log(initparam$pf)
    initparam$Cvec <- apply(initparam$data, 2, max) + 1
  }
  else if (scoretype == "usr") {
    if (is.null(usrpar$pctesttype)) {
      usrpar$pctesttype <- "usr"
    }
    initparam$pctesttype <- usrpar$pctesttype
    initparam <- usrscoreparameters(initparam, usrpar)
  }
  else if (scoretype == "mixed") {
    initparam$nbin <- mixedpar$nbin
    initparam$binpar <- scoreparameters("bde", data[, 1:mixedpar$nbin], 
                                        bdepar = bdepar, nodeslabels = nodeslabels[1:mixedpar$nbin], 
                                        weightvector = weightvector)
    initparam$gausspar <- scoreparameters("bge", data, bgnodes = c(1:mixedpar$nbin), 
                                          bgepar = bgepar, nodeslabels = nodeslabels, weightvector = weightvector)
  }
  attr(initparam, "class") <- "scoreparameters"
  return(initparam)
}
bidag_dagscore = function (scorepar, incidence) 
{
  if (scorepar$DBN) {
    stop("To calculate DBN score DBNscore should be used!")
  }
  n <- ncol(scorepar$data)
  if (scorepar$bgn == 0) {
    mainnodes <- c(1:scorepar$n)
  }
  else {
    mainnodes <- c(1:n)[-scorepar$bgnodes]
  }
  P_local <- numeric(n)
  for (j in mainnodes) {
    parentnodes <- which(incidence[, j] == 1)
    P_local[j] <- DAGcorescore(j, parentnodes, scorepar$n, 
                               scorepar)
  }
  return(sum(P_local))
}
bidag_dagcorescore = function (j, parentnodes, n, param) 
{
    if (param$DBN) {
        if (param$stationary) {
            internalparents <- parentnodes[which(parentnodes <= 
                param$nsmall)]
            corescore <- DAGcorescore(j, parentnodes, param$n + 
                param$nsmall, param$otherslices) + DAGcorescore(j, 
                internalparents, param$n, param$firstslice)
        }
        else {
            corescore <- 0
            for (i in 1:(length(param$paramsets) - 1)) {
                corescore <- corescore + DAGcorescore(j, parentnodes, 
                  param$n + param$nsmall, param$paramsets[[i]])
            }
            internalparents <- parentnodes[which(parentnodes <= 
                param$nsmall)]
            corescore <- corescore + DAGcorescore(j, internalparents, 
                param$n, param$paramsets[[length(param$paramsets)]])
        }
    }
    else if (param$MDAG) {
        corescore <- 0
        for (i in 1:length(param$paramsets)) {
            corescore <- corescore + DAGcorescore(j, parentnodes, 
                param$n, param$paramsets[[i]])
        }
    }
    else if (param$type == "bge") {
        TN <- param$TN
        awpN <- param$awpN
        scoreconstvec <- param$scoreconstvec
        lp <- length(parentnodes)
        awpNd2 <- (awpN - n + lp + 1)/2
        A <- TN[j, j]
        switch(as.character(lp), `0` = {
            corescore <- scoreconstvec[lp + 1] - awpNd2 * log(A)
        }, `1` = {
            D <- TN[parentnodes, parentnodes]
            logdetD <- log(D)
            B <- TN[j, parentnodes]
            logdetpart2 <- log(A - B^2/D)
            corescore <- scoreconstvec[lp + 1] - awpNd2 * logdetpart2 - 
                logdetD/2
            if (!is.null(param$logedgepmat)) {
                corescore <- corescore - param$logedgepmat[parentnodes, 
                  j]
            }
        }, `2` = {
            D <- TN[parentnodes, parentnodes]
            detD <- dettwobytwo(D)
            logdetD <- log(detD)
            B <- TN[j, parentnodes]
            logdetpart2 <- log(dettwobytwo(D - (B) %*% t(B)/A)) + 
                log(A) - logdetD
            corescore <- scoreconstvec[lp + 1] - awpNd2 * logdetpart2 - 
                logdetD/2
            if (!is.null(param$logedgepmat)) {
                corescore <- corescore - sum(param$logedgepmat[parentnodes, 
                  j])
            }
        }, {
            D <- as.matrix(TN[parentnodes, parentnodes])
            choltemp <- chol(D)
            logdetD <- 2 * log(prod(choltemp[(lp + 1) * c(0:(lp - 
                1)) + 1]))
            B <- TN[j, parentnodes]
            logdetpart2 <- log(A - sum(backsolve(choltemp, B, 
                transpose = TRUE)^2))
            corescore <- scoreconstvec[lp + 1] - awpNd2 * logdetpart2 - 
                logdetD/2
            if (!is.null(param$logedgepmat)) {
                corescore <- corescore - sum(param$logedgepmat[parentnodes, 
                  j])
            }
        })
    }
    else if (param$type == "bde") {
        lp <- length(parentnodes)
        noparams <- 2^lp
        chi <- param$chi
        scoreconstvec <- param$scoreconstvec
        switch(as.character(lp), `0` = {
            N1 <- sum(param$d1[, j])
            N0 <- sum(param$d0[, j])
            NT <- N0 + N1
            corescore <- scoreconstvec[lp + 1] + lgamma(N0 + 
                chi/(2 * noparams)) + lgamma(N1 + chi/(2 * noparams)) - 
                lgamma(NT + chi/noparams)
        }, `1` = {
            corescore <- scoreconstvec[lp + 1]
            summys <- param$data[, parentnodes]
            for (i in 1:noparams - 1) {
                totest <- which(summys == i)
                N1 <- sum(param$d1[totest, j])
                N0 <- sum(param$d0[totest, j])
                NT <- N0 + N1
                corescore <- corescore + lgamma(N0 + chi/(2 * 
                  noparams)) + lgamma(N1 + chi/(2 * noparams)) - 
                  lgamma(NT + chi/noparams)
            }
            if (!is.null(param$logedgepmat)) {
                corescore <- corescore - param$logedgepmat[parentnodes, 
                  j]
            }
        }, {
            summys <- colSums(2^(c(0:(lp - 1))) * t(param$data[, 
                parentnodes]))
            N1s <- collectC(summys, param$d1[, j], noparams)
            N0s <- collectC(summys, param$d0[, j], noparams)
            NTs <- N1s + N0s
            corescore <- scoreconstvec[lp + 1] + sum(lgamma(N0s + 
                chi/(2 * noparams))) + sum(lgamma(N1s + chi/(2 * 
                noparams))) - sum(lgamma(NTs + chi/noparams))
            if (!is.null(param$logedgepmat)) {
                corescore <- corescore - sum(param$logedgepmat[parentnodes, 
                  j])
            }
        })
    }
    else if (param$type == "bdecat") {
        lp <- length(parentnodes)
        chi <- param$chi
        corescore <- param$scoreconstvec[lp + 1]
        Cj <- param$Cvec[j]
        switch(as.character(lp), `0` = {
            Cp <- 1
            summys <- rep(0, nrow(param$data))
        }, `1` = {
            Cp <- param$Cvec[parentnodes]
            summys <- param$data[, parentnodes]
            if (!is.null(param$logedgepmat)) {
                corescore <- corescore - param$logedgepmat[parentnodes, 
                  j]
            }
        }, {
            Cp <- prod(param$Cvec[parentnodes])
            summys <- colSums(cumprod(c(1, param$Cvec[parentnodes[-lp]])) * 
                t(param$data[, parentnodes]))
            if (!is.null(param$logedgepmat)) {
                corescore <- corescore - sum(param$logedgepmat[parentnodes, 
                  j])
            }
        })
        if (!is.null(param$weightvector)) {
            Ns <- collectCcatwt(summys, param$data[, j], param$weightvector, 
                Cp, Cj)
        }
        else {
            Ns <- collectCcat(summys, param$data[, j], Cp, Cj)
        }
        NTs <- rowSums(Ns)
        corescore <- corescore + sum(lgamma(Ns + chi/(Cp * Cj))) - 
            sum(lgamma(NTs + chi/Cp)) + Cp * lgamma(chi/Cp) - 
            (Cp * Cj) * lgamma(chi/(Cp * Cj))
    }
    else if (param$type == "usr") {
        corescore <- usrDAGcorescore(j, parentnodes, n, param)
    }
    return(corescore)
}
bidag_usrscoreparameters = function (initparam, usrpar = list(pctesttype = "usrCItest", 
                                                          edgepf = 2, edgepmat = NULL, chi = 0.5, delta = NULL, eta = NULL)) 
{
  if (is.null(usrpar$chi)) {
    usrpar$chi <- 0.5
  }
  if (is.null(usrpar$edgepf)) {
    usrpar$edgepf <- 2
  }
  initparam$chi <- usrpar$chi
  initparam$pf <- usrpar$edgepf
  if (is.null(usrpar$delta)) {
    usrpar$delta <- 100 * initparam$chi
  }
  if (is.null(usrpar$eta)) {
    usrpar$eta <- 10 * initparam$chi
  }
  initparam$delta <- usrpar$delta
  initparam$eta <- usrpar$eta
  if (is.null(initparam$weightvector)) {
    initparam$N <- nrow(initparam$data)
    initparam$d1 <- initparam$data
    initparam$d0 <- (1 - initparam$data)
  }
  else {
    initparam$N <- sum(initparam$weightvector)
    initparam$d1 <- initparam$data * initparam$weightvector
    initparam$d0 <- (1 - initparam$data) * initparam$weightvector
  }
  maxparents <- ncol(initparam$data) - 1
  initparam$scoreconstvec <- rep(0, maxparents + 1)
  if (is.null(usrpar$edgepmat)) {
    initparam$logedgepmat <- NULL
  }
  else {
    initparam$logedgepmat <- log(usrpar$edgepmat)
  }
  initparam$scoreconstvec <- lgamma(initparam$chi/2) + lgamma((1 + 
                                                                 initparam$delta) * initparam$chi/4) - 3 * lgamma(initparam$chi/4) - 
    lgamma(initparam$delta * initparam$chi/4) - c(0:maxparents) * 
    log(initparam$pf)
  initparam$scoreconstvec[1] <- lgamma((1 + initparam$eta) * 
                                         initparam$chi/2) - lgamma(initparam$chi/2) - lgamma(initparam$eta * 
                                                                                               initparam$chi/2)
  initparam
}
bidag_usrDAGcorescore = function (j, parentnodes, n, param) 
{
  lp <- length(parentnodes)
  chi <- param$chi
  scoreconstvec <- param$scoreconstvec
  switch(as.character(lp), `0` = {
    N1 <- sum(param$d1[, j])
    N0 <- sum(param$d0[, j])
    NT <- N0 + N1
    corescore <- scoreconstvec[lp + 1] + lgamma(N0 + param$eta * 
                                                  chi/2) + lgamma(N1 + chi/2) - lgamma(NT + (1 + param$eta) * 
                                                                                         chi/2)
  }, `1` = {
    corescore <- scoreconstvec[lp + 1]
    summys <- param$data[, parentnodes]
    for (i in 0:1) {
      totest <- which(summys == i)
      N1 <- sum(param$d1[totest, j])
      N0 <- sum(param$d0[totest, j])
      NT <- N0 + N1
      if (i == 0) {
        corescore <- corescore + lgamma(N0 + param$delta * 
                                          chi/4) + lgamma(N1 + chi/4) - lgamma(NT + (1 + 
                                                                                       param$delta) * chi/4)
      } else {
        corescore <- corescore + lgamma(N0 + chi/4) + 
          lgamma(N1 + chi/4) - lgamma(NT + chi/2)
      }
    }
    if (!is.null(param$logedgepmat)) {
      corescore <- corescore - param$logedgepmat[parentnodes, 
                                                 j]
    }
  }, {
    summys <- 1 * (rowSums(param$data[, parentnodes]) == 
                     lp)
    N1s <- collectC(summys, param$d1[, j], 2)
    N0s <- collectC(summys, param$d0[, j], 2)
    NTs <- N1s + N0s
    corescore <- scoreconstvec[lp + 1] + sum(lgamma(N0s + 
                                                      c(param$delta, 1) * chi/4)) + sum(lgamma(N1s + chi/4)) - 
      sum(lgamma(NTs + c(1 + param$delta, 2) * chi/4))
    if (!is.null(param$logedgepmat)) {
      corescore <- corescore - sum(param$logedgepmat[parentnodes, 
                                                     j])
    }
  })
  corescore
}
bidag_iterativeMCMC <- function (scorepar, MAP = TRUE, posterior = 0.5, softlimit = 9, 
    hardlimit = 12, alpha = 0.05, gamma = 1, verbose = TRUE, 
    chainout = FALSE, scoreout = FALSE, cpdag = FALSE, mergetype = "skeleton", 
    iterations = NULL, moveprobs = NULL, stepsave = NULL, startorder = NULL, 
    accum = FALSE, compress = TRUE, plus1it = NULL, startspace = NULL, 
    blacklist = NULL, addspace = NULL, scoretable = NULL, alphainit = NULL) 
{
    if (is.null(moveprobs)) {
        prob1 <- 99
        if (scorepar$nsmall > 3) {
            prob1 <- round(6 * 99 * scorepar$nsmall/(scorepar$nsmall^2 + 
                10 * scorepar$nsmall - 24))
        }
        prob1 <- prob1/100
        moveprobs <- c(prob1, 0.99 - prob1, 0.01)
        moveprobs <- moveprobs/sum(moveprobs)
        moveprobs <- c(moveprobs[c(1, 2)], 0, moveprobs[3])
    }
    if (is.null(iterations)) {
        if (scorepar$nsmall < 26) {
            iterations <- 25000
        }
        else {
            iterations <- (3.5 * scorepar$nsmall * scorepar$nsmall * 
                log(scorepar$nsmall)) - (3.5 * scorepar$nsmall * 
                scorepar$nsmall * log(scorepar$nsmall))%%1000
        }
    }
    if (is.null(stepsave)) {
        stepsave <- floor(iterations/1000)
    }
    ordercheck <- checkstartorder(startorder, varnames = scorepar$labels.short, 
        mainnodes = scorepar$mainnodes, bgnodes = scorepar$static, 
        DBN = scorepar$DBN, split = scorepar$split)
    if (ordercheck$errorflag) {
        stop(ordercheck$message)
    }
    else {
        startorder <- ordercheck$order
    }
    if (scorepar$DBN) {
        if (!is.null(blacklist)) {
            blacklist <- DBNbacktransform(blacklist, scorepar)
        }
        if (!is.null(startspace)) {
            startspace <- DBNbacktransform(startspace, scorepar)
        }
        if (!is.null(addspace)) {
            addspace <- DBNbacktransform(addspace, scorepar)
        }
        if (scorepar$split) {
            if (scorepar$MDAG) {
                param1 <- scorepar$paramsets[[scorepar$nsets]]
                param2 <- scorepar$paramsets[[1]]
                param2$paramsets <- scorepar$paramsets[1:(scorepar$nsets - 
                  1)]
                param2$MDAG <- TRUE
            }
            else {
                param1 <- scorepar$firstslice
                param2 <- scorepar$otherslices
            }
            if (scoreout | !is.null(scoretable)) {
                cat("option scoreout always equals FALSE for DBNs with samestruct=FALSE, scoretable parameter is ignored \n")
            }
            cat("learning initial structure...\n")
            result.init <- iterativeMCMCplus1(param = param1, 
                iterations, stepsave, plus1it = plus1it, MAP = MAP, 
                posterior = posterior, alpha = alpha, cpdag = cpdag, 
                moveprobs = moveprobs, softlimit = softlimit, 
                hardlimit = hardlimit, startspace = startspace$init, 
                blacklist = blacklist$init, gamma = gamma, verbose = verbose, 
                chainout = chainout, scoreout = FALSE, mergecp = mergetype, 
                addspace = addspace$init, scoretable = NULL, 
                startorder = startorder$init, accum = accum, 
                alphainit = alphainit, compress = compress)
            cat("learning transition structure...\n")
            result.trans <- iterativeMCMCplus1(param = param2, 
                iterations, stepsave, plus1it = plus1it, MAP = MAP, 
                posterior = posterior, alpha = alpha, cpdag = cpdag, 
                moveprobs = moveprobs, softlimit = softlimit, 
                hardlimit = hardlimit, startspace = startspace$trans, 
                blacklist = blacklist$trans, gamma = gamma, verbose = verbose, 
                chainout = chainout, scoreout = FALSE, mergecp = mergetype, 
                addspace = addspace$trans, scoretable = NULL, 
                startorder = startorder$trans, accum = accum, 
                alphainit = alphainit, compress = compress)
            result <- mergeDBNres.it(result.init, result.trans, 
                scorepar)
        }
        else {
            result <- iterativeMCMCplus1(param = scorepar, iterations, 
                stepsave, plus1it = plus1it, MAP = MAP, posterior = posterior, 
                alpha = alpha, cpdag = cpdag, moveprobs = moveprobs, 
                softlimit = softlimit, hardlimit = hardlimit, 
                startspace = startspace, blacklist = blacklist, 
                gamma = gamma, verbose = verbose, chainout = chainout, 
                scoreout = scoreout, mergecp = mergetype, addspace = addspace, 
                scoretable = scoretable, startorder = startorder, 
                accum = accum, alphainit = alphainit, compress = compress)
        }
    }
    else {
        result <- iterativeMCMCplus1(param = scorepar, iterations, 
            stepsave, plus1it = plus1it, MAP = MAP, posterior = posterior, 
            alpha = alpha, cpdag = cpdag, moveprobs = moveprobs, 
            softlimit = softlimit, hardlimit = hardlimit, startspace = startspace, 
            blacklist = blacklist, gamma = gamma, verbose = verbose, 
            chainout = chainout, scoreout = scoreout, mergecp = mergetype, 
            addspace = addspace, scoretable = scoretable, startorder = startorder, 
            accum = accum, compress = compress)
    }
    result$info <- list()
    result$info$DBN <- scorepar$DBN
    if (scorepar$DBN) {
        result$info$nsmall <- scorepar$nsmall
        result$info$bgn <- scorepar$bgn
        result$info$split <- scorepar$split
    }
    result$info$algo <- "iterative order MCMC"
    if (is.null(startspace)) {
        result$info$spacealgo <- "PC"
    }
    else {
        result$info$spacealgo <- "user defined matrix"
    }
    result$info$iterations <- iterations
    result$info$plus1it <- length(result$max)
    result$info$samplesteps <- floor(iterations/stepsave) + 1
    if (MAP) {
        result$info$sampletype <- "MAP"
    }
    else {
        result$info$sampletype <- "sample"
        result$info$threshold <- posterior
    }
    result$info$fncall <- match.call()
    attr(result, "class") <- "iterativeMCMC"
    return(result)
}
bidag_iterativeMCMCplus1 <- function (param, iterations, stepsave, plus1it = NULL, MAP = TRUE, 
    posterior = 0.5, startorder = NULL, moveprobs, softlimit = 9, 
    hardlimit = 14, chainout = FALSE, scoreout = FALSE, startspace = NULL, 
    blacklist = NULL, gamma = 1, verbose = FALSE, alpha = NULL, 
    cpdag = FALSE, mergecp = "skeleton", addspace = NULL, scoretable = NULL, 
    accum, alphainit = NULL, compress = TRUE) 
{
    n <- param$n
    nsmall <- param$nsmall
    matsize <- ifelse(param$DBN, n + nsmall, n)
    objsizes <- list()
    maxlist <- list()
    maxobj <- list()
    updatenodeslist <- list()
    MCMCtraces <- list()
    if (!param$DBN) {
        if (param$bgn != 0) {
            updatenodes <- c(1:n)[-param$bgnodes]
        }
        else {
            updatenodes <- c(1:n)
        }
    }
    else {
        updatenodes <- c(1:nsmall)
    }
    if (is.null(blacklist)) {
        blacklist <- matrix(0, nrow = matsize, ncol = matsize)
    }
    diag(blacklist) <- 1
    if (!is.null(param$bgnodes)) {
        for (i in param$bgnodes) {
            blacklist[, i] <- 1
        }
    }
    if (!is.null(scoretable)) {
        startskel <- scoretable$adjacency
        blacklist <- scoretable$blacklist
        scoretable <- scoretable$tables
    }
    else {
        if (is.null(startspace)) {
            startspace <- definestartspace(alpha, param, cpdag = cpdag, 
                algo = "pc", alphainit = alphainit)
        }
        startskeleton <- 1 * (startspace & !blacklist)
        if (!is.null(addspace)) {
            startskel <- 1 * ((addspace | startskeleton) & !blacklist)
        }
        else {
            startskel <- startskeleton
        }
    }
    blacklistparents <- list()
    for (i in 1:matsize) {
        blacklistparents[[i]] <- which(blacklist[, i] == 1)
    }
    if (verbose) {
        cat(paste("maximum parent set size is", max(apply(startskel, 
            2, sum))), "\n")
    }
    if (max(apply(startskel, 2, sum)) > hardlimit) {
        stop("the size of maximal parent set is higher that the hardlimit; redifine the search space or increase the hardlimit!")
    }
    maxorder <- startorder
    ptab <- listpossibleparents.PC.aliases(startskel, isgraphNEL = FALSE, 
        n, updatenodes)
    if (verbose) {
        cat("core space defined, score table are being computed \n")
        flush.console()
    }
    parenttable <- ptab$parenttable
    aliases <- ptab$aliases
    numberofparentsvec <- ptab$numberofparentsvec
    numparents <- ptab$numparents
    plus1lists <- PLUS1(matsize, aliases, updatenodes, blacklistparents)
    rowmaps <- parentsmapping(parenttable, numberofparentsvec, 
        n, updatenodes)
    if (is.null(scoretable)) {
        scoretable <- scorepossibleparents.PLUS1(parenttable = parenttable, 
            plus1lists = plus1lists, n = n, param = param, updatenodes = updatenodes, 
            rowmaps, numparents, numberofparentsvec)
    }
    posetparenttable <- poset(parenttable, numberofparentsvec, 
        rowmaps, n, updatenodes)
    if (MAP == TRUE) {
        maxmatrices <- posetscoremax(posetparenttable, scoretable, 
            numberofparentsvec, rowmaps, n, plus1lists = plus1lists, 
            updatenodes)
    }
    else {
        bannedscore <- poset.scores(posetparenttable, scoretable, 
            ptab$numberofparentsvec, rowmaps, n, plus1lists = plus1lists, 
            ptab$numparents, updatenodes)
    }
    oldadj <- startskeleton
    i <- 1
    if (is.null(plus1it)) 
        plus1it <- 100
    while (length(updatenodes) > 0 & i <= plus1it) {
        if (i > 1) {
            newptab <- listpossibleparents.PC.aliases(newadj, 
                isgraphNEL = FALSE, n, updatenodes)
            parenttable[updatenodes] <- newptab$parenttable[updatenodes]
            aliases[updatenodes] <- newptab$aliases[updatenodes]
            numberofparentsvec[updatenodes] <- newptab$numberofparentsvec[updatenodes]
            numparents[updatenodes] <- newptab$numparents[updatenodes]
            newplus1lists <- PLUS1(matsize, aliases, updatenodes, 
                blacklistparents)
            plus1lists$mask[updatenodes] <- newplus1lists$mask[updatenodes]
            plus1lists$parents[updatenodes] <- newplus1lists$parents[updatenodes]
            plus1lists$aliases[updatenodes] <- newplus1lists$aliases[updatenodes]
            rowmaps[updatenodes] <- parentsmapping(parenttable, 
                numberofparentsvec, n, updatenodes)[updatenodes]
            scoretable[updatenodes] <- scorepossibleparents.PLUS1(parenttable, 
                plus1lists, n, param, updatenodes, rowmaps, numparents, 
                numberofparentsvec)[updatenodes]
            posetparenttable[updatenodes] <- poset(parenttable, 
                numberofparentsvec, rowmaps, n, updatenodes)[updatenodes]
            if (MAP) {
                newmaxmatrices <- posetscoremax(posetparenttable, 
                  scoretable, numberofparentsvec, rowmaps, n, 
                  plus1lists = plus1lists, updatenodes)
                maxmatrices$maxmatrix[updatenodes] <- newmaxmatrices$maxmatrix[updatenodes]
                maxmatrices$maxrow[updatenodes] <- newmaxmatrices$maxrow[updatenodes]
            }
            else {
                newbannedscore <- poset.scores(posetparenttable, 
                  scoretable, numberofparentsvec, rowmaps, n, 
                  plus1lists = plus1lists, numparents, updatenodes)
                bannedscore[updatenodes] <- newbannedscore[updatenodes]
            }
            if (verbose) {
                cat(paste("search space expansion", i, "\n"))
                flush.console()
            }
        }
        else {
            if (verbose) {
                cat(paste("score tables completed, iterative MCMC is running", 
                  "\n"))
                flush.console()
            }
        }
        if (MAP) {
            MCMCresult <- orderMCMCplus1max(n, nsmall, startorder, 
                iterations, stepsave, moveprobs, parenttable, 
                scoretable, aliases, numparents, rowmaps, plus1lists, 
                maxmatrices, numberofparentsvec, gamma = gamma, 
                bgnodes = param$bgnodes, matsize = matsize, chainout = chainout, 
                compress = compress)
        }
        else {
            MCMCresult <- orderMCMCplus1(n, nsmall, startorder, 
                iterations, stepsave, moveprobs, parenttable, 
                scoretable, aliases, numparents, rowmaps, plus1lists, 
                bannedscore, numberofparentsvec, gamma = gamma, 
                bgnodes = param$bgnodes, matsize = matsize, chainout = TRUE, 
                compress = compress)
        }
        MCMCtraces$DAGscores[[i]] <- MCMCresult$DAGscores
        if (chainout) {
            if (param$DBN) {
                MCMCtraces$incidence[[i]] <- lapply(MCMCresult$incidence, 
                  function(x) DBNtransform(x, param = param))
                MCMCtraces$orders[[i]] <- lapply(MCMCresult$orders, 
                  order2var, varnames = param$firstslice$labels)
            }
            else {
                MCMCtraces$incidence[[i]] <- lapply(MCMCresult$incidence, 
                  function(x) assignLabels(x, param$labels))
                MCMCtraces$orders[[i]] <- lapply(MCMCresult$orders, 
                  order2var, varnames = param$labels)
            }
            MCMCtraces$orderscores[[i]] <- MCMCresult$orderscores
        }
        maxobj <- storemaxMCMC(MCMCresult, param)
        maxlist[[i]] <- maxobj
        maxN <- which.max(MCMCresult$DAGscores)
        if (i > 1) {
            if (maxobj$score > maxscore) {
                maxDAG <- maxobj$DAG
                maxorder <- maxobj$order
                maxscore <- maxobj$score
                maxit <- i
            }
        }
        else {
            maxDAG <- maxobj$DAG
            maxscore <- maxobj$score
            maxorder <- maxobj$order
            maxit <- 1
        }
        if (MAP) {
            newadj <- newspacemap(n, startskeleton, oldadj, softlimit, 
                hardlimit, blacklist, maxdag = MCMCresult$maxdag, 
                mergetype = mergecp, accum = accum)
        }
        else {
            newadj <- newspaceskel(n, startskeleton, oldadj, 
                softlimit, hardlimit, posterior, blacklist, MCMCtrace = MCMCresult[[1]], 
                mergetype = mergecp)
        }
        updatenodes <- which(apply(newadj == oldadj, 2, all) == 
            FALSE)
        updatenodeslist[[i]] <- updatenodes
        if (is.null(plus1it)) {
            oldadj <- newadj
        }
        else if (i < plus1it) {
            oldadj <- newadj
        }
        else if (!scoreout) {
            oldadj <- newadj
        }
        startorder <- c(MCMCresult$orders[[maxN]], param$bgnodes)
        i <- i + 1
    }
    addedge <- sum(newadj) - sum(startskeleton)
    result <- list()
    if (scoreout) {
        if (chainout) {
            output <- 4
        }
        else {
            output <- 3
        }
    }
    else {
        if (chainout) {
            output <- 2
        }
        else {
            output <- 1
        }
    }
    result$maxtrace <- maxlist
    result$DAG <- maxobj$DAG
    result$CPDAG <- Matrix(graph2m(dag2cpdag(m2graph(result$DAG))), 
        sparse = TRUE)
    result$score <- maxobj$score
    result$maxorder <- maxobj$order
    result$trace <- MCMCtraces$DAGscores
    MCMCtraces$DAGscores <- NULL
    if (param$DBN) {
        result$startspace <- DBNtransform(startskeleton, param)
        result$endspace <- DBNtransform(oldadj, param)
    }
    else {
        result$startspace <- startskeleton
        result$endspace <- oldadj
    }
    switch(as.character(output), `1` = {
    }, `2` = {
        result$traceadd <- MCMCtraces
    }, `3` = {
        result$scoretable <- list()
        result$scoretable$adjacency <- result$endspace
        result$scoretable$tables <- scoretable
        result$scoretable$blacklist <- blacklist
        attr(result$scoretable, "class") <- "scorespace"
    }, `4` = {
        result$traceadd <- MCMCtraces
        result$scoretable <- list()
        result$scoretable$adjacency <- result$endspace
        result$scoretable$tables <- scoretable
        result$scoretable$blacklist <- blacklist
        attr(result$scoretable, "class") <- "scorespace"
    })
    return(result)
}
BiDAG_definestartspace <- function (alpha, param, cpdag = FALSE, algo = "pc", alphainit = NULL) 
{
    if (is.null(alphainit)) {
        alphainit <- alpha
    }
    local_type <- param$type
    if (local_type == "usr") {
        if (param$pctesttype %in% c("bde", "bge", "bdecat")) {
            local_type <- param$pctesttype
        }
    }
    if (param$DBN) {
        if (param$stationary) {
            othersliceskel <- definestartspace(alpha, param$otherslices, 
                cpdag = FALSE, algo = "pc")
            firstsliceskel <- definestartspace(alphainit, param$firstslice, 
                cpdag = FALSE, algo = "pc")
            startspace <- othersliceskel
            startspace[param$intstr$rows, param$intstr$cols] <- 1 * 
                (startspace[param$intstr$rows, param$intstr$cols] | 
                  firstsliceskel[param$intstr$rows, param$intstr$cols])
        }
        else {
            skels <- list()
            skels[[1]] <- definestartspace(alphainit, param$paramsets[[1]], 
                cpdag = FALSE, algo = "pc")
            startspace <- skels[[1]]
            for (i in 2:(length(param$paramsets) - 1)) {
                skels[[i]] <- definestartspace(alpha, param$paramsets[[i]], 
                  cpdag = FALSE, algo = "pc")
                startspace <- 1 * (skels[[i]] | startspace)
            }
            firstsliceskel <- definestartspace(alphainit, param$paramsets[[length(param$paramsets)]], 
                cpdag = FALSE, algo = "pc")
            startspace[param$intstr$rows, param$intstr$cols] <- 1 * 
                (startspace[param$intstr$rows, param$intstr$cols] | 
                  firstsliceskel[param$intstr$rows, param$intstr$cols])
        }
    }
    else {
        if (local_type == "bde") {
            if (cpdag) {
                pc.skel <- pc(suffStat = list(d1 = param$d1, 
                  d0 = param$d0, data = param$data), indepTest = weightedbinCItest, 
                  alpha = alpha, labels = colnames(param$data), 
                  verbose = FALSE)
            }
            else {
                pc.skel <- pcalg::skeleton(suffStat = list(d1 = param$d1, 
                  d0 = param$d0, data = param$data), indepTest = weightedbinCItest, 
                  alpha = alpha, labels = colnames(param$data), 
                  verbose = FALSE)
            }
        }
        else if (local_type == "bdecat") {
            if (cpdag) {
                pc.skel <- pc(suffStat = param, indepTest = weightedcatCItest, 
                  alpha = alpha, labels = colnames(param$data), 
                  verbose = FALSE)
            }
            else {
                pc.skel <- pcalg::skeleton(suffStat = param, 
                  indepTest = weightedcatCItest, alpha = alpha, 
                  labels = colnames(param$data), verbose = FALSE)
            }
        }
        else if (local_type == "bge") {
            if (is.null(param$weightvector)) {
                cormat <- cor(param$data)
                N <- nrow(param$data)
            }
            else {
                N <- sum(param$weightvector)
                cormat <- cov.wt(param$data, wt = param$weightvector, 
                  cor = TRUE)$cor
            }
            if (cpdag) {
                pc.skel <- pcalg::pc(suffStat = list(C = cormat, 
                  n = N), indepTest = gaussCItest, alpha = alpha, 
                  labels = colnames(param$data), skel.method = "stable", 
                  verbose = FALSE)
            }
            else {
                pc.skel <- pcalg::skeleton(suffStat = list(C = cormat, 
                  n = N), indepTest = gaussCItest, alpha = alpha, 
                  labels = colnames(param$data), method = "stable", 
                  verbose = FALSE)
            }
        }
        else if (local_type == "usr") {
            if (cpdag) {
                pc.skel <- pc(suffStat = param, indepTest = usrCItest, 
                  alpha = alpha, labels = colnames(param$data), 
                  verbose = FALSE)
            }
            else {
                pc.skel <- pcalg::skeleton(suffStat = param, 
                  indepTest = usrCItest, alpha = alpha, labels = colnames(param$data), 
                  verbose = FALSE)
            }
        }
        g <- pc.skel@graph
        startspace <- 1 * (graph2m(g))
    }
    return(startspace)
}