#bidag fns
bidag_score_param = function (scoretype = c("bge", "bde", "bdecat", "usr"), data, 
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
bidag_score_dag = function (scorepar, incidence) 
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
bidag_usr_scoreparam = function (initparam, usrpar = list(pctesttype = "usrCItest", 
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
bidag_usr_score_dag = function (j, parentnodes, n, param) 
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
####################################################################
# dp bidag
####################################################################
bidag_bge = function (scoretype = c("bge", "bde", "bdecat", "usr"), data, 
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