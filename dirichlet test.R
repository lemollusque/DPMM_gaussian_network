library(dirichletprocess)

its <- 500

# Gaussian
y <- rt(200, 3) + 2 #generate sample data 
dp <- DirichletProcessGaussian(y) 
dp <- Fit(dp, its)



faithfulTransformed <- scale(faithful$waiting) 
dp <- DirichletProcessGaussian(faithfulTransformed) 
dp <- Fit(dp, its) 
plot(dp, likelihood =TRUE, single=TRUE) 
plot(dp, data_method="hist")

# f <- LikelihoodFunction(dp)(dp$data)
pf = PosteriorFunction(dp)(dp$data)
# log likelihood conditional on DAG
loglikelihood = sum(log(pf))

plot(dp$data, f)
plot(dp$data, pf)



# By default, the Fit function updates the cluster allocation, the cluster 
# parameters and then the α parameter. 
# In some rare cases, updating α every iteration can delay convergence. 
# Instead, the user can instead choose to update α every 10 iterations.
dp <- DirichletProcessGaussian(y) 
samples <- list() 
for(s in seq_len(1000)){ 
  dp <- ClusterComponentUpdate(dp) 
  dp <- ClusterParameterUpdate(dp) 
  
  if(s %% 10 == 0) { 
    dp <- UpdateAlpha(dp) 
  } 
  samples[[s]] <- list() 
  samples[[s]]$phi <- dp$clusterParameters 
  samples[[s]]$weights <- dp$weights 
}

# Multivariate normal
faithfulTrans <- scale(faithful)
dp <- DirichletProcessMvnormal(faithfulTrans)
dp <- Fit(dp, 200)
plot(dp) 

md = dp$mixingDistribution
new_md = MixingDistribution(distribution="mvnormal", 
                   priorParameters = md$priorParameters, 
                   conjugate="conjugate")
dp <- DirichletProcessCreate(faithfulTrans, new_md) 
dp <- Initialise(dp)
dp <- Fit(dp, 200) 
pf <- PosteriorFrame(dp, faithfulTrans, 200) 




y <- rnorm(10)
dp <- DirichletProcessGaussian(y)
dp <- Fit(dp, 5)
f <- LikelihoodFunction(dp)
PosteriorFunction(dp)
plot(-2:2, f(-2:2))
plot(dp)
