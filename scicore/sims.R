packages <- c(
  "BiDAG",
  "dirichletprocess",
  "dplyr",
  "ggplot2",
  "foreach",
  "doFuture",
  "future",
  "progressr",
  "doRNG",
  "mvtnorm"
)

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

source("comparison_algs.R")
source("dualPC.R")
source("dao.R")
source("fns.R")
insertSource("fns.R", package = "BiDAG")

init.seed <- 100
iter <- 100
dual <- FALSE

# dirichlet params
alpha_prior <- c(2, 4)
dp_iter <- 500
dp_fits <- 4
burnin <- 400
L <- 20

param_grid <- expand.grid(
  N = c(100),
  n = 4,
  d = c(0, 1, 2, 5, 10),
  bge.par = 1
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
    "Sims/",
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
    .packages = c("BiDAG", "matrixStats", "dirichletprocess", "dplyr", "mclust", "mvtnorm")
  ) %dorng% {
    
    source("comparison_algs.R")
    source("dualPC.R")
    source("dao.R")
    source("fns.R")
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
    job_seed <- make_job_seed(init.seed, k)
    set.seed(job_seed)
    
    g <- er_dag(n)
    g <- sf_out(g)
    truegraph <- randomize_graph(g)
    
    model1 <- corr(truegraph)
    model2 <- corr(truegraph)
    
    X1 <- simulate(model1$B, model1$O, N / 2)
    X2 <- simulate(model2$B, model2$O, N / 2)
    
    v <- rnorm(ncol(X2))
    v <- v / sqrt(sum(v^2))
    shift <- d * v
    X2 <- sweep(X2, 2, shift, "+")
    
    data <- standardize(rbind(X1, X2))
    if (is.null(colnames(data))) {
      colnames(data) <- paste0("v", seq_len(ncol(data)))
    }
    
    # prepare dirichlet gamma list
    Gamma_list <- list()
    for (f in seq_len(dp_fits)) {
      dp <- DirichletProcessMvnormal(data, alphaPriors = alpha_prior)
      dp <- Fit(dp, dp_iter, progressBar = FALSE)
      
      Gamma_sample <- dp_membership_probs(dp, burnin, L)
      Gamma_list <- add_membershipp(
        Gamma_list,
        Gamma_sample,
        child = colnames(data)[1],
        parents = colnames(data)[2:ncol(data)],
        active = TRUE
      )
    }
    
    dp_usrpar <- list(
      pctesttype = "bge",
      am = bge.par,
      membershipp_list = Gamma_list,
      edgepf = 1
    )
    
    # search spaces
    DP.searchspace <- set.searchspace(data, dual, "DP", usrpar = dp_usrpar)
    bge.searchspace <- set.searchspace(data, dual, "bge", bge.par)
    
    iter_results <- data.frame()
    
    bge.fit <- bge.partition.mcmc(bge.searchspace, order = FALSE)
    iter_results <- compare_results(
      bge.fit, c(bge.par, "BGe, partition"), iter_results, truegraph
    )
    
    bge.fit <- bge.partition.mcmc(bge.searchspace, order = TRUE)
    iter_results <- compare_results(
      bge.fit, c(bge.par, "BGe, order"), iter_results, truegraph
    )
    
    dp.fit <- DP.partition.mcmc(DP.searchspace, order = FALSE)
    iter_results <- compare_results(
      dp.fit, c(bge.par, "DP, partition"), iter_results, truegraph
    )
    
    dp.fit <- DP.partition.mcmc(DP.searchspace, order = TRUE)
    iter_results <- compare_results(
      dp.fit, c(bge.par, "DP, order"), iter_results, truegraph
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

