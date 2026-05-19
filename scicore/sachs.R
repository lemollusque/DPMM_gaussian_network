library(BiDAG)
library(matrixStats)
library(dirichletprocess)
library(dplyr)
library(ggplot2)
library(foreach)
library(doFuture)
library(future)
library(parallelly)
library(mclust)
library(progressr)
library(doRNG)
library(mvtnorm)
library(readxl)

source("comparison_algs.R")
source("dualPC.R")
source("dao.R")
source("fns.R")
insertSource("fns.R", package = "BiDAG")

sachs.data <- read.csv("Sachs/2005_sachs_2_cd3cd28icam2_log_std.csv")
sachs.data <- as.matrix(sachs.data)
trueDAGbn <- read.csv("Sachs/sachs.csv")
trueDAGbn <- as.matrix(trueDAGbn)

bge.par = 0.01
dual <- TRUE 
# dirichlet params
dp_iter <- 100
burnin <- 50
L <- 10
dp_fits <- 1

# dp settings
dp_usrpar <- list(
  pctesttype = "bge",
  am = bge.par,
  dp_iter = dp_iter,
  dp_burnin = burnin,
  dp_n_sample = L,
  dp_fits = dp_fits
)

init.seed <- 100
iter <- 20

dir.create("Sachs/parallel_results", showWarnings = FALSE, recursive = TRUE)

make_file_name <- function(i) {
  paste0("Sachs/parallel_results/sachs_rep_", sprintf("%03d", i), ".rds")
}

make_job_seed <- function(init.seed, i) {
  init.seed + i
}

n_cores <- max(1, availableCores() - 1)

plan(multisession, workers = n_cores)
registerDoFuture()
registerDoRNG(init.seed)

handlers(global = TRUE)
handlers("txtprogressbar")

with_progress({
  p <- progressor(steps = iter)
  
  foreach(
    i = seq_len(iter),
    .packages = c(
      "BiDAG", "matrixStats", "dirichletprocess", "dplyr",
      "mclust", "mvtnorm"
    )
  ) %dopar% {
    
    source("comparison_algs.R")
    source("dualPC.R")
    source("dao.R")
    source("fns.R")
    insertSource("fns.R", package = "BiDAG")
    
    file_name <- make_file_name(i)
    
    if (file.exists(file_name)) {
      p(sprintf("skip %d", i))
      return(NULL)
    }
    
    set.seed(make_job_seed(init.seed, i))
    
    iter_results <- data.frame()
    
    DP.searchspace <- set.searchspace(
      sachs.data,
      dual,
      "DP",
      usrpar = dp_usrpar
    )
    
    bge.searchspace <- set.searchspace(
      sachs.data,
      dual,
      "bge",
      bge.par
    )
    
    # DP score, partition
    dp.fit <- DP.partition.mcmc(
      DP.searchspace,
      order = FALSE,
      iterations = 1200
    )
    dp.edgep <- post.edges(dp.fit)
    iter_results <- compare_results(
      dp.fit,
      c(dp.edgep, "DP, partition"),
      iter_results,
      trueDAG
    )
    
    # DP score, order
    dp.fit <- DP.partition.mcmc(
      DP.searchspace,
      order = TRUE,
      iterations = 1200
    )
    dp.edgep <- post.edges(dp.fit)
    iter_results <- compare_results(
      dp.fit,
      c(dp.edgep, "DP, order"),
      iter_results,
      trueDAG
    )
    
    # BGe score, partition
    bge.fit <- bge.partition.mcmc(
      bge.searchspace,
      order = FALSE,
      iterations = 1200
    )
    bge.edgep <- post.edges(bge.fit)
    iter_results <- compare_results(
      bge.fit,
      c(bge.edgep, "BGe, partition"),
      iter_results,
      trueDAG
    )
    
    # BGe score, order
    bge.fit <- bge.partition.mcmc(
      bge.searchspace,
      order = TRUE,
      iterations = 1200
    )
    bge.edgep <- post.edges(bge.fit)
    iter_results <- compare_results(
      bge.fit,
      c(bge.edgep, "BGe, order"),
      iter_results,
      trueDAG
    )
    
    colnames(iter_results) <- c(
      "ESHD", "eTP", "eFP", "TPR", "FPR_P", "time",
      "ErktoAkt", "ErktoPKA", "Scorefn", "graph"
    )
    
    iter_results$rep <- i
    iter_results$job_seed <- make_job_seed(init.seed, i)
    
    saveRDS(iter_results, file_name)
    
    p(sprintf("done %d", i))
  }
})