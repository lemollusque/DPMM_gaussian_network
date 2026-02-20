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
dual <- F  # use dualPC
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
############################### compare all dags ##############################
# List all DAGs with n nodes
all.dags <- list()
adj <- matrix(0, nrow = n, ncol = n,
              dimnames=list(vars, vars))
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

# Compute true posterior (and save all scores)
true.post <- rep(NA, dag.counter)

for(d in 1:dag.counter) {
  print(d)
  dag <- all.dags[[d]]
  true.post[d] = BiDAG::DAGscore(usr_score_param, dag)
}
# Order true posterior
true.order <- order(true.post, decreasing = T)
true.post <- true.post[true.order]
true.p <- exp(true.post - logSumExp(true.post))
all.dags <- all.dags[true.order]


#compute shd for every graph
shd <- rep(NA, dag.counter)
count_edges <- rep(NA, dag.counter)
for(k in 1:dag.counter) {
  print(k)
  dag <- all.dags[[k]]
  camparison = compareDAGs(dag, truegraph)
  shd[k] = camparison["SHD"]
  count_edges[k] = sum(dag)
}
shd.order =  order(shd, decreasing = F) 
a <- true.post[shd.order]
b <- true.p[shd.order]


par(mfrow = c(2,1))
plot(a, type="l", ylab="posterior")
plot(b, type="l", ylab="softmax")




top5_idx <- order(true.post, decreasing = T)[1:5]
all.dags[top5_idx]
true.post[top5_idx]
