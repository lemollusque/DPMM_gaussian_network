# preprocessing etc of the Sachs data
source("./isachssetup.R")

# load libraries
library(BiDAG)
library(Bestie)

library(matrixStats)
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
library(BayesFactor)
library(matrixStats)
library(BNPmix)

# Use BiDAG with intervention scoring
insertSource("./usrscorefns.R", package = "BiDAG")
# Use Bestie with intervention scoring
insertSource("./usrparamfns.R", package = "Bestie")

# load causal pipeline taken and adapted from https://github.com/annlia/causalpipe
source("itoyDAGfunctionsSachs.R")
source("intfns.R")
source("dualPC.R")


library(data.table) # for last
library(DiagrammeR) # for making DAG plot
library(DiagrammeRsvg) ## for exporting svg for plotting to file
library(rsvg) ## for converting svg to png

sachs.data <- scale(data)
bgn <- ncol(Imat)

init.seed <- 234
set.seed(init.seed)

bge.par = 0.01
# dirichlet params
dp_usrpar <- list(
  pctesttype = "bge",
  am = bge.par,
  dp_prior = list(strength = 0.05, discount = 0),
  dp_mcmc = list(niter = 4000, nburn = 3999, model="LS"),
  dp_n_sample = 100,
  dp_fits = 1,
  dp_fitspace = "full"
)
output <- list(out_param = TRUE, out_type = "FULL")  


cor_mat <- cor(sachs.data)
fitspace <- dual_pc(cor_mat, nrow(sachs.data), alpha = 0.05, skeleton = T)
child = "p38"
parents <- names(which(fitspace[ , child] == 1))
children <- names(which(fitspace[child, ] == 1))
neighbour <- setdiff(unique(c(parents, children)), child)

dp_data = sachs.data[,c(child, neighbour)]
fit_em <- Mclust(
  data = dp_data,
  G = 6
)
n
hard_clusters <- fit_em$classification
table(hard_clusters)

clusters = unique(hard_clusters)

nDAGs <- 50
nSeeds <- 50
batch <- 100 + 1:nSeeds
labels4plot <- colnames(sachs.data) 
nNodes <- length(labels4plot)

plan(multisession, workers = min(length(batch), availableCores() - 1))
registerDoFuture()

for (cluster in clusters){
  data = sachs.data[hard_clusters == cluster,, drop = FALSE]
  ImatCluster = Imat[hard_clusters == cluster,, drop = FALSE]
  exps <- mgcv::uniquecombs(ImatCluster)
  n_exps <- nrow(exps)
  if (is.null(n_exps)) {
    n_exps <- 1L
  }
  scoreObject_ibge <- scoreparameters(scoretype = "usr", data = data, usrpar = list(pctesttype = "bge", Imat = ImatCluster, am = bge.par))
  dname = paste0("cluster", cluster)
  foreach(
    seednumber = batch,
    .packages = c("BiDAG", "Bestie", "data.table", "mvtnorm")
  ) %dorng% {
    
    # load causal pipeline taken and adapted from https://github.com/annlia/causalpipe
    source("itoyDAGfunctionsSachs.R")
    source("intfns.R")
    
    print(paste("Seed is", seednumber))
    
    load(file = paste0("./saveout_subpopulation/dagdraw", nNodes, "seed", seednumber, "", ".RData"))
    if (n_exps == 1L) {
      sampledDAGs <- lapply(sampledDAGs, function(A) {
        A[-(1:bgn), -(1:bgn), drop = FALSE]
      })
      
    }
    causalMats <- DAGintervention(sampledDAGs, scoreObject_ibge, sample=TRUE)
    save(causalMats,
         file=paste0("./saveout_subpopulation/effects", nNodes, "seed", seednumber, dname, ".RData"))
    
    TRUE
  }
}


# p38 -> jnk analysis
p38 <- which(colnames(sachs.data) == "p38") 
jnk <- which(colnames(sachs.data) == "JNK") 

eff_df <- data.frame()

for (cl in clusters) {
  
  dname <- paste0("cluster", cl)
  ImatCluster = Imat[hard_clusters == cl,, drop = FALSE]
  exps <- mgcv::uniquecombs(ImatCluster)
  n_exps <- nrow(exps)
  if (is.null(n_exps)) {
    print(1)
    n_exps <- 1L
  }
  
  alleffs <- vector("list", nDAGs * length(batch))
  for (nlevel in seq_along(batch)) {
    seednumber <- batch[nlevel]
    
    load(file = paste0(
      "./saveout_subpopulation/effects",
      nNodes, "seed", seednumber, dname, ".RData"
    ))
    
    alleffs[1:nDAGs + (nlevel - 1) * nDAGs] <- causalMats
  }
  if (n_exps == 1){
    eff <- sapply(alleffs, function(x) x[p38, jnk])
  }
  else{
    eff <- sapply(alleffs, function(x) x[p38 + bgn, jnk + bgn])
  }
  
  eff <- eff[eff != 0]
  
  eff_df <- rbind(
    eff_df,
    data.frame(
      effect = eff,
      cluster = paste0("cluster ", cl)
    )
  )
}

# cluster contribution by number of effect samples
cluster_weights <- as.data.frame(
  prop.table(table(hard_clusters))
)

names(cluster_weights) <- c("cluster", "weight")
cluster_weights$cluster <- paste0(
  "cluster ",
  cluster_weights$cluster
)

eff_df <- eff_df %>%
  left_join(cluster_weights, by = "cluster")

ggplot(eff_df, aes(x = effect, colour = cluster, fill = cluster)) +
  geom_density(
    aes(weight = weight / n),
    alpha = 0.25,
    linewidth = 0.8
  ) +
  geom_density(
    data = eff_df,
    aes(x = effect),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 1
  ) +
  coord_cartesian(ylim = c(0, 7)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(
    x = "P38 → Jnk",
    y = "Density"
  ) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(vjust = 0.5, hjust = 1),
    panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  ) +
  guides(
    colour = "none",
    fill = guide_legend(title = NULL)
  )


# compare densities
# DP-BGe
dname <- ""
alleffs <- vector("list", nDAGs * length(batch))
for (nlevel in seq_along(batch)) {
  seednumber <- batch[nlevel]
  
  load(file = paste0(
    "./saveout_subpopulation/effects",
    nNodes, "seed", seednumber, dname, ".RData"
  ))
  
  alleffs[1:nDAGs + (nlevel - 1) * nDAGs] <- causalMats
}

eff <- sapply(alleffs, function(x) x[p38, jnk])
alleff_dp <- eff[eff != 0]

# BGe
dname = "bge"
scoreObject_ibge <- scoreparameters(scoretype = "usr", data = sachs.data, usrpar = list(pctesttype = "bge", Imat = Imat, am = bge.par))
foreach(
  seednumber = batch,
  .packages = c("BiDAG", "Bestie", "data.table", "mvtnorm")
) %dorng% {
  
  source("itoyDAGfunctionsSachs.R")
  source("intfns.R")
  
  print(paste("Seed is", seednumber))
  
  load(file = paste0("./saveout_subpopulation/dagdraw", nNodes, "seed", seednumber, "", ".RData"))
  causalMats <- DAGintervention(sampledDAGs, scoreObject_ibge, sample=TRUE)
  save(causalMats,
       file=paste0("./saveout_subpopulation/effects", nNodes, "seed", seednumber, dname, ".RData"))
  
  TRUE
}
alleffs <- vector("list", nDAGs * length(batch))
for (nlevel in seq_along(batch)) {
  seednumber <- batch[nlevel]
  
  load(file = paste0(
    "./saveout_subpopulation/effects",
    nNodes, "seed", seednumber, dname, ".RData"
  ))
  
  alleffs[1:nDAGs + (nlevel - 1) * nDAGs] <- causalMats
}
eff <- sapply(alleffs, function(x) x[p38+bgn, jnk+bgn])
alleff_bge <- eff[eff != 0]

method_cols <- c(
  "iBGe" = "#F8766D",
  "DP-iBGe" = "#619CFF",
  "Cluster mixture" = "black"
)

plot_df <- rbind(
  data.frame(effect = alleff_dp, source = "DP-iBGe"),
  data.frame(effect = alleff_bge, source = "iBGe"),
  data.frame(effect = eff_df$effect, source = "Cluster mixture")
)
plot_df$source <- factor(
  plot_df$source,
  levels = c("iBGe", "DP-iBGe", "Cluster mixture")
)

ggplot(plot_df,
       aes(x = effect,
           colour = source,
           linetype = source)) +
  geom_density(linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_colour_manual(values = method_cols) +
  scale_linetype_manual(
    values = c(
      "iBGe" = "solid",
      "DP-iBGe" = "solid",
      "Cluster mixture" = "solid"
    )
  ) +
  coord_cartesian(ylim = c(0, 7)) +
  theme_bw() +
  labs(
    x = "P38 → Jnk",
    y = "Density",
    colour = NULL,
    linetype = NULL
  ) +
  theme(
    legend.position = "bottom",
    panel.grid.major.x = element_blank()
  )
