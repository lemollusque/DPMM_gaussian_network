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

library(data.table) # for last
library(DiagrammeR) # for making DAG plot
library(DiagrammeRsvg) ## for exporting svg for plotting to file
library(rsvg) ## for converting svg to png

source("toyDAGfunctionsSachs.R")

source("comparison_algs.R")
source("Fourier_fns.R")
source("dualPC.R")
source("dao.R")
source("fns.R")
insertSource("fns.R", package = "BiDAG")

sachs.data <- read.csv("Sachs/2005_sachs_2_cd3cd28icam2_log_std.csv")
sachs.data <- as.matrix(sachs.data)

init.seed <- 234
set.seed(init.seed)

bge.par = 0.01
# dirichlet params
dp_usrpar <- list(
  pctesttype = "bge",
  am = bge.par,
  dp_prior = list(strength = 1, discount = 0),
  dp_mcmc = list(niter = 4000, nburn = 3000, model="LS"),
  dp_n_sample = 100,
  dp_fits = 1,
  dp_fitspace = "full"
)
output <- list(out_param = TRUE, out_type = "FULL")  


cor_mat <- cor(sachs.data)
fitspace <- dual_pc(cor_mat, nrow(sachs.data), alpha = 0.05, skeleton = T)
child = "Jnk"
parents <- names(which(fitspace[ , child] == 1))
children <- names(which(fitspace[child, ] == 1))
neighbour <- setdiff(unique(c(parents, children)), child)

dp_data = sachs.data[,c(child, neighbour)]
fit <- PYdensity(y = dp_data, mcmc = dp_usrpar$dp_mcmc, prior = dp_usrpar$dp_prior, output = output)

part <- partition(fit)
hard_clusters <- part$partitions[1, ]
table(hard_clusters)

data = sachs.data[hard_clusters == 1,]
dname = "cluster1"

nDAGs <- 50
nSeeds <- 50
batch <- 100 + 1:nSeeds
labels4plot <- colnames(sachs.data) 
nNodes <- length(labels4plot)

plan(multisession, workers = min(length(batch), availableCores() - 1))
registerDoFuture()

foreach(
  seednumber = batch,
  .packages = c("BiDAG", "Bestie", "data.table", "mvtnorm")
) %dorng% {
  
  source("fns.R")
  insertSource("fns.R", package = "BiDAG")
  source("toyDAGfunctionsSachs.R")
  
  print(paste("Seed is", seednumber))
  
  
  load(file = paste0("./Sachs/saveout_subpopulation/dagdraw", nNodes, "seed", seednumber, "", ".RData"))
  
  bge.searchspace <- set.searchspace(data, "bge", bge.par)
  causalMats <- DAGintervention(sampledDAGs, bge.searchspace$score, sample=TRUE)
  
  save(causalMats,
       file=paste0("./Sachs/saveout_subpopulation/effects", nNodes, "seed", seednumber, dname, ".RData"))
  
  TRUE
}

nSeeds <- length(batch)
alldigraphs <- vector("list", nDAGs * nSeeds) # to store the graphs
alleffs <- vector("list", nDAGs * nSeeds) # to store the matrices of effects
for (nlevel in 1:nSeeds) {
  seednumber <- batch[nlevel]
  ## Retrieve sampled DAGs - DAG chain
  load(file = paste0("./Sachs/saveout_subpopulation/dagdraw", nNodes, "seed", seednumber, "", ".RData"))
  alldigraphs[1:nDAGs + (nlevel - 1) * nDAGs] <- sampledDAGs # remove the starting point
  ## Retrieve estimated effects
  load(file = paste0("./Sachs/saveout_subpopulation/effects", nNodes, "seed", seednumber, dname, ".RData"))
  alleffs[1:nDAGs + (nlevel - 1) * nDAGs] <- causalMats # remove the starting point
}
data4plot = list(alldigraphs=alldigraphs, alleffs=alleffs)

pdf(paste0("Sachs/SachsEffects_", dname, ".pdf"), width = 6, height = 6)
plotEffects(effects4plot = data4plot$alleffs, xmargs = c(0.1, 0.3), label_size = 1.5,
            sortlabs = 1:11, title_text = "")
dev.off()




# jnk -> p38 analysis
jnk <- which(colnames(sachs.data) == "Jnk")
p38 <- which(colnames(sachs.data) == "P38")

par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
for (cl in 1:6) {
  
  dname <- paste0("cluster", cl)
  
  nSeeds <- length(batch)
  alleffs <- vector("list", nDAGs * nSeeds)
  
  for (nlevel in 1:nSeeds) {
    seednumber <- batch[nlevel]
    
    load(file = paste0(
      "./Sachs/saveout_subpopulation/effects",
      nNodes, "seed", seednumber, dname, ".RData"
    ))
    
    alleffs[1:nDAGs + (nlevel - 1) * nDAGs] <- causalMats
  }
  
  eff <- sapply(alleffs, function(x) x[jnk, p38])
  eff <- eff[eff != 0]
  
  hist(
    eff,
    breaks = 50,
    main = dname,
    xlab = "Jnk → P38 effect",
    col = "grey80",
    border = "white"
  )
}

cluster = hard_clusters

eff_df <- data.frame()

for (cl in 1:6) {
  
  dname <- paste0("cluster", cl)
  alleffs <- vector("list", nDAGs * length(batch))
  
  for (nlevel in seq_along(batch)) {
    seednumber <- batch[nlevel]
    
    load(file = paste0(
      "./Sachs/saveout_subpopulation/effects",
      nNodes, "seed", seednumber, dname, ".RData"
    ))
    
    alleffs[1:nDAGs + (nlevel - 1) * nDAGs] <- causalMats
  }
  
  eff <- sapply(alleffs, function(x) x[jnk, p38])
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
cluster_weights <- eff_df %>%
  group_by(cluster) %>%
  summarise(n = n()) %>%
  mutate(weight = n / sum(n))

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
    x = "Jnk → P38",
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
alleff_dp <- sapply(data4plot$alleffs, function(x) x[jnk, p38])
alleff_dp <- alleff_dp[alleff_dp != 0]

# BGe
alleff_bge <- sapply(data4plot_bge$alleffs, function(x) x[jnk, p38])
alleff_bge <- alleff_bge[alleff_bge != 0]

method_cols <- c(
  "BGe" = "#F8766D",
  "DP" = "#619CFF",
  "Cluster mixture" = "black"
)

plot_df <- rbind(
  data.frame(effect = alleff_dp, source = "DP"),
  data.frame(effect = alleff_bge, source = "BGe"),
  data.frame(effect = eff_df$effect, source = "Cluster mixture")
)
plot_df$source <- factor(
  plot_df$source,
  levels = c("BGe", "DP", "Cluster mixture")
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
      "BGe" = "solid",
      "DP" = "solid",
      "Cluster mixture" = "solid"
    )
  ) +
  coord_cartesian(ylim = c(0, 7)) +
  theme_bw() +
  labs(
    x = "Jnk → P38",
    y = "Density",
    colour = NULL,
    linetype = NULL
  ) +
  theme(
    legend.position = "bottom",
    panel.grid.major.x = element_blank()
  )
