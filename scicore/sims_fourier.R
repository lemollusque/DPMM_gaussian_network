packages <- c(
  "BiDAG",
  "dplyr",
  "ggplot2",
  "foreach",
  "doFuture",
  "future",
  "progressr",
  "doRNG",
  "mvtnorm",
  "BayesFactor",
  "matrixStats",
  "BNPmix"
)

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

source("comparison_algs.R")
source("dualPC.R")
source("Fourier_fns.R")
source("dao.R")
source("fns.R")
insertSource("fns.R", package = "BiDAG")

init.seed <- 234
iter <- 100
dp_fitspace = "full"

# dirichlet params
dp_iter <- 5000
burnin <- 3000
L <- 100

param_grid <- expand.grid(
  N = c(100, 200, 500, 1000),
  n = 10,
  d = c(0, 0.5, 1),
  bge.par = 0.01
)

sim_grid <- expand.grid(
  j = seq_len(nrow(param_grid)),
  i = seq_len(iter)
)

make_job_seed <- function(init.seed, k) {
  init.seed + k
}

make_file_name <- function(N, n, d, bge.par, i) {
  paste0(
    "Sims_fourier/",
    "N", N,
    "_n", n,
    "_d", d,
    "_bge", bge.par,
    "_rep", sprintf("%03d", i),
    ".rds"
  )
}

n_cores <- max(1, availableCores())

plan(multisession, workers = n_cores)
registerDoFuture()
registerDoRNG(init.seed)

handlers(global = TRUE)
handlers("txtprogressbar")

with_progress({
  p <- progressor(steps = nrow(sim_grid))
  
  foreach(
    k = seq_len(nrow(sim_grid)),
    .packages = c("BiDAG", "matrixStats", "dplyr", "mclust", "mvtnorm")
  ) %dopar% {
    
    source("comparison_algs.R")
    source("dualPC.R")
    source("dao.R")
    source("fns.R")
    source("Fourier_fns.R")
    insertSource("fns.R", package = "BiDAG")
    
    # set params
    j <- sim_grid$j[k]
    i <- sim_grid$i[k]
    
    N <- param_grid$N[j]
    n <- param_grid$n[j]
    d <- param_grid$d[j]
    bge.par <- param_grid$bge.par[j]
    
    file_name <- make_file_name(
      N = N,
      n = n,
      d = d,
      bge.par = bge.par,
      i = i
    )
    
    # resume logic: skip completed jobs
    if (file.exists(file_name)) {
      p(sprintf("skip %d", k))
      return(NULL)
    }
    
    # deterministic per-job seed
    job_seed <- make_job_seed(init.seed, i)
    set.seed(job_seed)
    
    myDAG <- pcalg::randomDAG(n, prob = 0.2, lB = 1, uB = 2) 
    trueDAG <- as(myDAG, "matrix")
    truegraph <- 1*(trueDAG != 0)
    data <- Fou_nldata(truegraph, N, lambda = d, noise.sd = 1, standardize = T)
    
    if (is.null(colnames(data))) {
      colnames(data) <- paste0("v", seq_len(ncol(data)))
    }
    
    
    iter_results <- data.frame()
    
    # save bge results
    bge.searchspace <- set.searchspace(data, "bge", bge.par)
    bge.fit <- bge.partition.mcmc(bge.searchspace, order = FALSE)
    iter_results <- compare_results(
      bge.fit, c(bge.par, "BGe, partition"), iter_results, detectable_truegraph
    )
    bge.fit <- bge.partition.mcmc(bge.searchspace, order = TRUE)
    iter_results <- compare_results(
      bge.fit, c(bge.par, "BGe, order"), iter_results, detectable_truegraph
    )
    
    # save dp results
    # dp settings
    dp_usrpar <- list(
      pctesttype = "bge",
      am = bge.par,
      dp_prior = list(strength = 1, discount = 0, model="LS"),
      dp_mcmc = list(niter = dp_iter, nburn = burnin),
      dp_n_sample = L,
      dp_fits = dp_fits,
      dp_fitspace = "full"
    )

    DP.searchspace <- set.searchspace(data, "DP", usrpar = dp_usrpar)
    dp.fit <- DP.partition.mcmc(DP.searchspace, order = FALSE)
    iter_results <- compare_results(
      dp.fit, c(bge.par, "DP, partition"), iter_results,  detectable_truegraph
    )
    dp.fit <- DP.partition.mcmc(DP.searchspace, order = TRUE)
    iter_results <- compare_results(
      dp.fit, c(bge.par, "DP, order"), iter_results,  detectable_truegraph
    )
    
    # save dp sub results
    # dp settings
    dp_usrpar <- list(
      pctesttype = "bge",
      am = bge.par,
      dp_prior = list(strength = 1, discount = 0, model="LS"),
      dp_mcmc = list(niter = dp_iter, nburn = burnin),
      dp_n_sample = L,
      dp_fits = dp_fits,
      dp_fitspace = "dual"
    )
    
    DP.searchspace <- set.searchspace(data, "DP", usrpar = dp_usrpar)
    dp.fit <- DP.partition.mcmc(DP.searchspace, order = FALSE)
    iter_results <- compare_results(
      dp.fit, c(bge.par, "DP dual, partition"), iter_results,  detectable_truegraph
    )
    dp.fit <- DP.partition.mcmc(DP.searchspace, order = TRUE)
    iter_results <- compare_results(
      dp.fit, c(bge.par, "DP dual, order"), iter_results,  detectable_truegraph
    )
    
    iter_results$N <- N
    iter_results$n <- n
    iter_results$d <- d
    iter_results$rep <- i
    iter_results$bge.par <- bge.par
    iter_results$job_id <- k
    iter_results$job_seed <- job_seed
    
    colnames(iter_results) <- c(
      "ESHD", "eTP", "eFP", "TPR", "FPR_P",
      "time", "parameter", "method", "graph",
      "N", "n", "d", "rep", "bge.par", "job_id", "job_seed"
    )
    
    saveRDS(iter_results, file_name)
    
    p(sprintf("done %d", k))
  }
})
