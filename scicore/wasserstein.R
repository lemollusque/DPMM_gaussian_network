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
  "BNPmix",
  "transport"
)

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

source("toyDAGfunctionsSachs.R")

source("comparison_algs.R")
source("Fourier_fns.R")
source("dualPC.R")
source("dao.R")
source("fns.R")
insertSource("fns.R", package = "BiDAG")

init.seed <- 234
iter <- 100
nDAGs <- 500

# dirichlet params
dp_model = "LS"
dp_strength = 0.1
dp_fits <- 1
dp_iter <- 5000
burnin <- 3000
L <- 20

param_grid <- expand.grid(
  N = c(100, 200, 500, 1000),
  n = 10,
  d = c(0, 1, 2, 5, 10),
  bge.par = 0.01
)

sim_grid <- expand.grid(
  j = seq_len(nrow(param_grid)),
  i = seq_len(iter)
)

make_file_name <- function(N, n, d, bge.par, i) {
  paste0(
    "wasserstein/",
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
    .packages = c(
      "BiDAG", "Bestie", "matrixStats", "dirichletprocess",
      "dplyr", "mclust", "mvtnorm", "BNPmix"
    )
  ) %dopar% {
    
    source("toyDAGfunctionsSachs.R")
    source("comparison_algs.R")
    source("Fourier_fns.R")
    source("dualPC.R")
    source("dao.R")
    source("fns.R")
    insertSource("fns.R", package = "BiDAG")
    
    j <- sim_grid$j[k]
    i <- sim_grid$i[k]
    
    N <- param_grid$N[j]
    n <- param_grid$n[j]
    d <- param_grid$d[j]
    bge.par <- param_grid$bge.par[j]
    
    file_name <- make_file_name(N, n, d, bge.par, i)
    
    if (file.exists(file_name)) {
      p(sprintf("skip %d", k))
      return(NULL)
    }
    
    job_seed <- init.seed + k
    set.seed(job_seed)
    
    myDAG <- pcalg::randomDAG(n, prob = 0.2, lB = 1, uB = 2)
    trueDAG <- as(myDAG, "matrix")
    truegraph <- 1 * (trueDAG != 0)
    
    sim <- simulate_bimodal(
      dag = t(truegraph),
      n = N,
      bimodal_sep = d,
      return_model = TRUE
    )
    
    truegraph <- t(sim$detectable_truegraph)
    
    data <- sim$data
    if (is.null(colnames(data))) {
      colnames(data) <- paste0("v", seq_len(ncol(data)))
    }
    
    B1 <- sim$model1$B
    B2 <- sim$model2$B
    
    Eff1 <- t(solve(diag(ncol(B1)) - B1))
    Eff2 <- t(solve(diag(ncol(B2)) - B2))
    
    rescale_effects <- function(Eff, sds) {
      D <- diag(sds)
      D %*% Eff %*% solve(D)
    }
    
    Eff1 <- rescale_effects(Eff1, sim$raw_sd)
    Eff2 <- rescale_effects(Eff2, sim$raw_sd)
    
    trueEffectSamples <- sampleTrueEffects(
      Eff1 = Eff1,
      Eff2 = Eff2,
      n1 = sim$n1,
      n2 = sim$n2,
      labels = colnames(truegraph),
      nSamples = 1000
    )
    
    tol <- 1e-10
    effectGraph <- 1 * ((abs(Eff1) > tol) | (abs(Eff2) > tol))
    diag(effectGraph) <- 0
    
    trueSamples <- makeTrueDAGSamples(
      truegraph = truegraph,
      nDAGs = nDAGs
    )
    
    ##################################################
    # DP Wasserstein
    ##################################################
    
    dp_usrpar <- list(
      pctesttype = "bge",
      am = bge.par,
      dp_prior = list(strength = dp_strength, discount = 0),
      dp_mcmc = list(niter = dp_iter, nburn = burnin, model = dp_model),
      dp_n_sample = L,
      dp_fits = dp_fits,
      dp_fitspace = "full"
    )
    
    DP.searchspace <- set.searchspace(
      data,
      "DP",
      usrpar = dp_usrpar
    )
    
    dp_alleffs <- computeEffects_mem(
      sampledDAGs = trueSamples,
      scoreObject = DP.searchspace$score,
      DP = TRUE
    )
    
    dp_res <- wasserstein_effect_avg(
      effects4plot = dp_alleffs,
      trueEffects = trueEffectSamples,
      effectGraph = effectGraph,
      sortlabs = seq_len(n)
    )
    
    # dual 
    dp_usrpar <- list(
      pctesttype = "bge",
      am = bge.par,
      dp_prior = list(strength = dp_strength, discount = 0),
      dp_mcmc = list(niter = dp_iter, nburn = burnin, model = dp_model),
      dp_n_sample = L,
      dp_fits = dp_fits,
      dp_fitspace = "dual"
    )
    
    DP.searchspace <- set.searchspace(
      data,
      "DP",
      usrpar = dp_usrpar
    )
    
    dp_dual_alleffs <- computeEffects_mem(
      sampledDAGs = trueSamples,
      scoreObject = DP.searchspace$score,
      DP = TRUE
    )
    
    dp_dual_res <- wasserstein_effect_avg(
      effects4plot = dp_dual_alleffs,
      trueEffects = trueEffectSamples,
      effectGraph = effectGraph,
      sortlabs = seq_len(n)
    )
    ##################################################
    # BGe Wasserstein
    ##################################################
    
    bge.searchspace <- set.searchspace(data, "bge", bge.par)
    
    bge_alleffs <- computeEffects_mem(
      sampledDAGs = trueSamples,
      scoreObject = bge.searchspace$score,
      DP = FALSE
    )
    
    bge_res <- wasserstein_effect_avg(
      effects4plot = bge_alleffs,
      trueEffects = trueEffectSamples,
      effectGraph = effectGraph,
      sortlabs = seq_len(n)
    )
    
    iter_results <- data.frame(
      method = c("DP", "DP dual", "BGe"),
      avgW = c(dp_res$avgW, dp_dual_res$avgW, bge_res$avgW),
      N = N,
      n = n,
      d = d,
      bge.par = bge.par,
      rep = i,
      job_id = k,
      job_seed = job_seed
    )
    iter_results <- rbind(
      data.frame(
        metric = "avgW",
        value = dp_res$avgW,
        method = "DP"
      ),
      data.frame(
        metric = "avgW",
        value = dp_dual_res$avgW,
        method = "DP dual"
      ),
      data.frame(
        metric = "avgW",
        value = bge_res$avgW,
        method = "BGe"
      )
    )
    
    iter_results$N <- N
    iter_results$n <- n
    iter_results$d <- d
    iter_results$rep <- i
    iter_results$bge.par <- bge.par
    iter_results$job_id <- k
    iter_results$job_seed <- job_seed
    
    saveRDS(iter_results, file_name)
    
    p(sprintf("done %d", k))
  }
})

files <- list.files("wasserstein", pattern = "\\.rds$", full.names = TRUE)
results <- bind_rows(lapply(files, readRDS))


ggplot(results, aes(x = method, y = value, color = method)) +
  geom_boxplot(aes(group = method), width = 0.6, outlier.shape = NA, linewidth = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.7, size = 0.5) +
  labs(x = NULL, y = "E-W") +
  facet_grid(
    N ~ d,
    scales = "free_y",
    labeller = labeller(
      N = function(x) paste("N =", x),
      d = function(x) paste("d =", x)
    )
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )


