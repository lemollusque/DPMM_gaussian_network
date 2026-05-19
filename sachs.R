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
library(GGally)
library(igraph)

source("comparison_algs.R")
source("dualPC.R")
source("dao.R")
source("fns.R")
insertSource("fns.R", package = "BiDAG")

sachs.data <- read.csv("Sachs/2005_sachs_2_cd3cd28icam2_log_std.csv")
sachs.data <- as.matrix(sachs.data)

# plot data
my_smooth <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_point(color = "black", alpha = 0.5, size = 0.6) +
    geom_smooth(
      method = "lm",
      se = FALSE,
      color = "red",
      linewidth = 0.8
    )
}

ggpairs(
  sachs.data,
  upper = list(continuous = wrap("cor", size = 4)),
  lower = list(continuous = my_smooth),
  diag = list(continuous = wrap("barDiag", bins = 30, alpha = 0.5))
) +
  theme_minimal()

trueDAG <- read.csv("Sachs/sachs.csv")
trueDAG <- as.matrix(trueDAG)

# plot DAG
g <- graph_from_adjacency_matrix(
  trueDAG,
  mode = "undirected",
  diag = FALSE
)

plot(
  g,
  layout = layout_with_fr(g),   # Fruchterman-Reingold
  vertex.size = 25,
  vertex.color = "white",
  vertex.frame.color = "black",
  vertex.label.cex = 0.8,
  edge.width = 1.5
)

set.seed(100)

results <- data.frame()

bge.par = 0.01
dual <- TRUE 
# dirichlet params
dp_iter <- 1000
burnin <- 500
L <- 100
dp_fits <- 4

# dp settings
dp_usrpar <- list(
  pctesttype = "bge",
  am = bge.par,
  dp_iter = dp_iter,
  dp_burnin = burnin,
  dp_n_sample = L,
  dp_fits = dp_fits
)

# search spaces
# DP.searchspace <- set.searchspace(sachs.data, dual, "DP", usrpar = dp_usrpar)
# saveRDS(DP.searchspace, "Sachs/dpsearchspace.rds")
DP.searchspace <- readRDS("Sachs/dpsearchspace.rds")

bge.searchspace = set.searchspace(sachs.data, dual, "bge", bge.par)

for (i in 1:100){
  print(i)
  # DP score, partition
  dp.fit <- DP.partition.mcmc(DP.searchspace, order = FALSE, iterations = 1200)
  dp.edgep <- post.edges(dp.fit)
  results <- compare_results(dp.fit, c(dp.edgep, "DP, partition"), results, trueDAGbn)
  
  # DP score, order
  dp.fit <- DP.partition.mcmc(DP.searchspace, order = TRUE, iterations = 1200)
  dp.edgep <- post.edges(dp.fit)
  results <- compare_results(dp.fit, c(dp.edgep, "DP, order"), results, trueDAGbn)
  
  # BGe score, partition
  bge.fit <-  bge.partition.mcmc(bge.searchspace, order = F, iterations = 1200)
  bge.edgep <- post.edges(bge.fit)
  results <- compare_results(bge.fit, c(bge.edgep, "BGe, partition"), results, trueDAGbn)
  
  # BGe score, order
  bge.fit <-  bge.partition.mcmc(bge.searchspace, order = T, iterations = 1200)
  bge.edgep <- post.edges(bge.fit)
  results <- compare_results(bge.fit, c(bge.edgep, "BGe, order"), results, trueDAGbn)
}

colnames(results) <- c("ESHD", "eTP", "eFP", "TPR", "FPR_P", "time",
                       "ErktoAkt", "ErktoPKA", "Scorefn", "graph")
saveRDS(results, "Sachs/Sachs_results.rds")

results <- readRDS("Sachs/Sachs_results.rds")

# plots resutts
keep_methods <- c("DP, partition", "DP, order", "BGe, partition", "BGe, order")

results_small <- results %>%
  filter(Scorefn %in% keep_methods) 

results_small <- results_small %>%
  mutate(ESHD = as.numeric(ESHD))
# add text medians
medians <- results_small %>%
  group_by(Scorefn, graph) %>%
  summarise(median_ESHD = median(ESHD), .groups = "drop")

ggplot(results_small, aes(x = Scorefn, y = ESHD, color = Scorefn)) +
  geom_boxplot(aes(group = Scorefn), width = 0.6, outlier.shape = NA, linewidth = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.7, size = 1.2) +
  geom_text(
    data = medians,
    aes(x = Scorefn, y = median_ESHD, label = round(median_ESHD,2)),
    color = "black",
    vjust = -0.7,
    size = 3
  ) +
  
  labs(x = NULL, y = "E-SHD") +
  theme_bw() +
  facet_grid( ~ graph, scales = "free_y") +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major.x = element_blank()
  )