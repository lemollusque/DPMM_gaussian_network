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
library(ggplotify)
library(patchwork)

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

nDAGs <- 50
nSeeds <- 50
batch <- 100 + 1:nSeeds
labels4plot <- colnames(sachs.data) 
nNodes <- length(labels4plot)

data4plot <- loadsamples(seeds=batch, nn=nNodes)

graph2plot <- dagviz(data4plot$alldigraphs, style_mat = matrix(1, 11, 11), title_text = "")
rsvg_png(charToRaw(export_svg(graph2plot)), "Sachs/SachsDAGs.png", width = 4000)


desired_order <- c(
  "Raf", "Mek", "Plcg", "PIP2", "PIP3",
  "Erk", "Akt", "PKA", "PKC", "P38", "Jnk"
)
sortlabs <- match(desired_order, colnames(data4plot$alleffs[[1]]))
pdf("Sachs/SachsEffects.pdf", width = 6, height = 6)
plotEffects(effects4plot = data4plot$alleffs, xmargs = c(0.1, 0.3), label_size = 1.5,
            sortlabs = sortlabs, title_text = "")
dev.off()

# subplots (make 4 subplots)
p1 <- as.ggplot(~plotEffects(effects4plot = data4plot$alleffs, xmargs = c(0.1, 0.3), label_size = 3, efflabs_size=2,
            sortlabs = sortlabs[1:2], title_text = ""))
p2 <- as.ggplot(~plotEffects(effects4plot = data4plot$alleffs, xmargs = c(0.1, 0.3), label_size = 3, efflabs_size=2,
                             sortlabs = sortlabs[3:5], title_text = ""))
p3 <- as.ggplot(~plotEffects(effects4plot = data4plot$alleffs, xmargs = c(0.1, 0.3), label_size = 3, efflabs_size=2,
                             sortlabs = sortlabs[6:8], title_text = ""))
p4 <- as.ggplot(~plotEffects(effects4plot = data4plot$alleffs, xmargs = c(0.1, 0.3), label_size = 3, efflabs_size=2,
                             sortlabs = sortlabs[9:11], title_text = ""))

style_effect_plot <- function(p) {
  p +
    theme(
      plot.background = element_rect(
        colour = "grey40",
        fill = "white",
        linewidth = 0.6
      ),
      plot.margin = margin(6, 6, 6, 6)
    )
}

p1 <- style_effect_plot(p1)
p2 <- style_effect_plot(p2)
p3 <- style_effect_plot(p3)
p4 <- style_effect_plot(p4)

combined_plot <-
  p1 + plot_spacer() +
  p2 + plot_spacer() +
  p3 + plot_spacer() +
  p4 +
  plot_layout(
    nrow = 1,
    widths = c(3, 0.15, 4, 0.15, 4, 0.15, 4)
  )

ggsave(
  "Sachs/SachsEffects_4_block.pdf",
  combined_plot,
  width = 16,
  height = 4
)
