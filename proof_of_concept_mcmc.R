library(dirichletprocess)
library(gRbase)
library(BiDAG)
library(matrixStats)
library(cowplot)
library(pcalg)
library(rstan)
library(bridgesampling)

source("fns.R")
source("Fourier_fns.R")
source("BayesStanFns.R")
source("sampling_fns.R")

set.seed(101)
lambda <- 1
n <- 4  # number of nodes
N <- 100  # number of samples

myDAG <- pcalg::randomDAG(n, prob = 0.5, lB = 1, uB = 2) 
trueDAG <- as(myDAG, "matrix")
truegraph <- 1*(trueDAG != 0)
data <- Fou_nldata(truegraph, N, lambda = lambda, noise.sd = 1, standardize = T) 

vars  <- c("x1","x2","x3","x4")
colnames(data) = vars

#----------------------  overwrite functions ----------------------------------
# replace BIDAG functions
unlockBinding("usrscoreparameters", asNamespace("BiDAG"))
assign("usrscoreparameters", usrscoreparameters, envir = asNamespace("BiDAG"))
lockBinding("usrscoreparameters", asNamespace("BiDAG"))

unlockBinding("usrDAGcorescore", asNamespace("BiDAG"))
assign("usrDAGcorescore", usrDAGcorescore, envir = asNamespace("BiDAG"))
lockBinding("usrDAGcorescore", asNamespace("BiDAG"))
#----------------------  end overwrite functions ------------------------------
# Initiate params for DP and BGe
n <- ncol(data)
alpha_mu <- 1          
alpha_w  <- n + alpha_mu + 1      
t <- alpha_mu * (alpha_w - n - 1) / (alpha_mu + 1)

# DP
g0Priors <- list(
  mu0    = rep(0, n),
  Lambda = diag(n) / t,   # T = (1/t) I
  kappa0 = alpha_mu,
  nu     = alpha_w
)

scaled_data = scale(data) 
n_iter = 100
burnin = 30
L = 10 # sample to take

Gamma_list <- list()
vars  <- c("x1","x2","x3","x4")
for (child in vars){
  parents <- vars[vars != child]
  dp_data = scaled_data[,c(child, parents)]
  dp <-  DirichletProcessMvnormal(dp_data, g0Priors)
  dp <- Fit(dp, n_iter)
  
  Gamma_sample <- dp_membership_probs(dp, n_iter, burnin, L)
  Gamma_list <- add_membershipp(Gamma_list, Gamma_sample, child=child, parents=parents)
}

# scoring
usr_score_param <- BiDAG::scoreparameters(scoretype = "usr", 
                                          data = scaled_data, 
                                          usrpar = list(pctesttype = "bge",
                                                        membershipp_list = Gamma_list,
                                                        am = alpha_mu, 
                                                        aw = alpha_w, 
                                                        T0scale = t,
                                                        edgepf = 1
                                          )
)

#----------------------------------- posterior -------------------------------
# search space
start <- Sys.time()
it_mcmc <- BiDAG::iterativeMCMC(scorepar = usr_score_param, 
                                    hardlimit = 14, 
                                    verbose = F, 
                                    scoreout = TRUE)
time <- Sys.time() - start

searchspace = list(score = usr_score_param, scoretable = it_mcmc$scoretable, DAG = it_mcmc$DAG, 
     maxorder = it_mcmc$maxorder, endspace = it_mcmc$endspace, time = time)

# partition mcmc
dp.mcmc <- function(searchspace, alpha = 0.05, 
                               order = FALSE, burnin = 0.33, iterations = 600) {
  start <- Sys.time()
  score <- searchspace$score
  
  if(order) {
    dp_mcmc_fit <- BiDAG::orderMCMC(score, MAP = FALSE, chainout = TRUE, alpha = alpha, 
                                startorder = searchspace$maxorder, scoretable = searchspace$scoretable,
                                startspace = searchspace$endspace, iterations = iterations, stepsave = 4)
  }
  else {
    dp_mcmc_fit <- BiDAG::partitionMCMC(score, alpha = alpha, startDAG = searchspace$DAG, 
                                    scoretable = searchspace$scoretable, startspace = searchspace$endspace,
                                    iterations = iterations, stepsave = 4)
  }
  toburn <- round(burnin * dp_mcmc_fit$info$samplesteps)
  
  dp_mcmc_fit$traceadd$incidence <- dp_mcmc_fit$traceadd$incidence[-(1:toburn)]
  time <- Sys.time() - start + searchspace$time
  dp_mcmc_fit$time <- as.numeric(time, units = "secs")
  
  return(dp_mcmc_fit)
}

dp_mcmc_fit = dp.mcmc(searchspace)
sampled_dags <- dp_mcmc_fit$traceadd$incidence

# Find index in all.dags of sampled dags 
# List all DAGs with n nodes
all.dags <- list()
adj <- matrix(0, nrow = n, ncol = n)
dag.counter <- 0
all.comb <- rep(list(c(0,1)), n*(n-1))
all.comb <- expand.grid(all.comb)  # all combinations outside of diagonal of adjacency matrix

for(i in 1:nrow(all.comb)) {
  adj[col(adj)!=row(adj)] <- as.numeric(all.comb[i, ])
  
  if(is.DAG(adj)) {
    dag.counter <- dag.counter + 1
    all.dags[[dag.counter]] <- adj
  }
}
all.vecdags <- lapply(all.dags, c)
post.indexes <- sapply(sampled_dags, function(x) which(all.vecdags %in% list(as.numeric(x))))
#BGe
dp_post.tab <- table(post.indexes)
dp_post.post <- as.numeric(dp_post.tab)/sum(dp_post.tab)
ind.dp_post <- as.numeric(names(dp_post.tab))

results <- data.frame(kl = KL_div(bge.post, ind.bge, true.p), 
                                     iter = iters[i], group = "DP")