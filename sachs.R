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

trueDAG <- read.csv("Sachs/sachs.csv")
trueDAG <- as.matrix(trueDAG)

set.seed(100)

bge.par = 0.01
dual <- TRUE 

results <- data.frame()

files <- list.files(
  "Sachs/parallel_searchspaces",
  pattern = "\\.rds$",
  full.names = TRUE
)

plan(multisession, workers = min(length(files), availableCores() - 1))
registerDoFuture()

results_list <- foreach(
  file = files,
  .packages = c("BiDAG", "Bestie", "data.table", "mvtnorm")
) %dorng% {
  
  source("fns.R")
  insertSource("fns.R", package = "BiDAG")
  source("toyDAGfunctionsSachs.R")
  
  cat("Running file:", file, "\n")
  
  DP.searchspace <- readRDS(file)
  bge.searchspace <- set.searchspace(sachs.data, dual, "bge", bge.par)
  
  res <- data.frame()
  
  dp.fit <- DP.partition.mcmc(DP.searchspace, order = FALSE, iterations = 1200)
  dp.edgep <- post.edges(dp.fit, weights = dp.fit$weights)
  res <- compare_results(dp.fit, c(dp.edgep, "DP, partition"), res, trueDAG)
  
  dp.fit <- DP.partition.mcmc(DP.searchspace, order = TRUE, iterations = 1200)
  dp.edgep <- post.edges(dp.fit, weights = dp.fit$weights)
  res <- compare_results(dp.fit, c(dp.edgep, "DP, order"), res, trueDAG)
  
  bge.fit <- bge.partition.mcmc(bge.searchspace, order = FALSE, iterations = 1200)
  bge.edgep <- post.edges(bge.fit)
  res <- compare_results(bge.fit, c(bge.edgep, "BGe, partition"), res, trueDAG)
  
  bge.fit <- bge.partition.mcmc(bge.searchspace, order = TRUE, iterations = 1200)
  bge.edgep <- post.edges(bge.fit)
  res <- compare_results(bge.fit, c(bge.edgep, "BGe, order"), res, trueDAG)
  
  res$file <- basename(file)
  res
}

results <- dplyr::bind_rows(results_list)

colnames(results) <- c("ESHD", "eTP", "eFP", "TPR", "FPR_P", "time",
                       "ErktoAkt", "ErktoPKA", "Scorefn", "graph", "file")

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