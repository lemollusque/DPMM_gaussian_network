setwd("C:/Users/BASTIAN/Desktop/code/git_DPMM_gaussian_network/DPMM_gaussian_network")

library(dirichletprocess)
library(gRbase)
library(BiDAG)
library(matrixStats)
library(questionr)
library(ggpubr)
library(cowplot)
library(pcalg)
library(rstan)
library(bridgesampling)

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



get_logl = function(obs, multinomial=T, iteration=500){
  obs = scale(obs)
  if(multinomial){
    dp = DirichletProcessMvnormal(obs)
  }
  else{
    dp = DirichletProcessGaussian(obs)
  }
  dp = Fit(dp, iteration)
  
  densities = Predictive(dp$mixingDistribution, obs)
  logd = log(densities)
  logl = sum(logd)
  
  return(logl)
}


##############################################################
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

# Compute true posterior (and save all scores)
true.post <- rep(NA, dag.counter)
parent.scores <- data.frame(sets = character(0), newscore = character(0))  
parent.scores <- rep(list(parent.scores), n)

for(k in 1:dag.counter) {
  print(k)
  dag <- all.dags[[k]]
  curr_score <- 0
  
  for(x in 1:n) {
    print(x)
    set <- parent.scores[[x]]$sets
    pax <- dag[,x]
    pax.str <- paste(pax, collapse = "")  # parents of x
    check <- is.na(match(set, pax.str))  # check if parent set is already there
    
    if(all(check)) {  # new parent set
      
      if(sum(pax) == 0) {  # Compute local score
        print("no pax")
        loc_score <- get_logl(data[,x], multinomial=F) 
      }
      else if(sum(pax) == 1){
        print("1 pax")
        logl_pax = get_logl(data[ ,which(pax==1)], multinomial=F)
        logl_xpax = get_logl(data[ ,c(x, which(pax==1))])
        loc_score = logl_xpax - logl_pax
      }
      else {
        print("more pax")
        logl_pax = get_logl(data[ ,which(pax==1)])
        logl_xpax = get_logl(data[ ,c(x, which(pax==1))])
        loc_score = logl_xpax - logl_pax
      }
      parent.scores[[x]][length(set)+1, ] <- c(pax.str, loc_score)
    }
    
    else {  # fetch score from parent.scores
      ind <- which(check == F)  # index of pax in set
      loc_score <- as.numeric(parent.scores[[x]]$newscore[ind])
    }
    curr_score <- curr_score + loc_score  # build score
  }
  true.post[k] <- curr_score
}

true.p <- exp(true.post - logSumExp(true.post))


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

plot(a, type="l")
plot(b, type="l")
