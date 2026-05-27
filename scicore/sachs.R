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

bge.par = 0.01
start_type = "full"
# dirichlet params
dp_iter <- 200
burnin <- 100
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

init.seed <- 234
iter <- 20

dir.create("Sachs/parallel_searchspaces", showWarnings = FALSE, recursive = TRUE)


make_dp_searchspace_file <- function(i) {
  paste0("Sachs/parallel_searchspaces/DP_searchspace_rep_", sprintf("%03d", i), ".rds")
}

make_job_seed <- function(init.seed, i) {
  init.seed + i
}

n_cores <- max(1, availableCores())

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
    
    dp_searchspace_file <- make_dp_searchspace_file(i)
    
    if (file.exists(dp_searchspace_file)) {
      p(sprintf("skip %d", i))
      return(NULL)
    }
    
    set.seed(make_job_seed(init.seed, i))
    
    iter_results <- data.frame()
    
    DP.searchspace <- set.searchspace(
      sachs.data,
      start_type,
      "DP",
      usrpar = dp_usrpar
    )
    saveRDS(DP.searchspace, dp_searchspace_file)
    
    
    p(sprintf("done %d", i))
  }
})