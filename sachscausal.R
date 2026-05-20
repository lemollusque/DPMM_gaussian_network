# load libraries
library(BiDAG)
library(Bestie)
source("fns.R")
insertSource("fns.R", package = "BiDAG")

# Use BiDAG with intervention scoring
# insertSource("../usrscorefns.R", package = "BiDAG")
# # Use Bestie with intervention scoring
# insertSource("../usrparamfns.R", package = "Bestie")

# load causal pipeline taken and adapted from https://github.com/annlia/causalpipe
source("toyDAGfunctionsSachs.R")

library(data.table) # for last
library(DiagrammeR) # for making DAG plot
library(DiagrammeRsvg) ## for exporting svg for plotting to file
library(rsvg) ## for converting svg to png

sachs.data <- read.csv("Sachs/2005_sachs_2_cd3cd28icam2_log_std.csv")
sachs.data <- as.matrix(sachs.data)

trueDAG <- read.csv("Sachs/sachs.csv")
trueDAG <- as.matrix(trueDAG)

nDAGs <- 50
nSeeds <- 1
batch <- 100 + 1:nSeeds
labels4plot <- colnames(sachs.data) 
nNodes <- length(labels4plot)



for(seednumber in batch){
  timing <- proc.time()
  print(paste("Seed is", seednumber))
  DP.searchspace <- readRDS("Sachs/dpsearchspace_sachs_2.rds")
  sampleDAGs(inData=sachs.data,
             scoreObject = DP.searchspace$score,
             nDigraphs=nDAGs, 
             seed=seednumber)
  computeEffects(n=nNodes, seed=seednumber)
  print(proc.time() - timing)
}

data4plot <- loadsamples(seeds=batch, nn=nNodes)

graph2plot <- dagviz(data4plot$alldigraphs, style_mat = matrix(1, 11, 11), title_text = "")
rsvg_png(charToRaw(export_svg(graph2plot)), "Sachs/SachsDAGs.png", width = 4000)

pdf("SachsEffects.pdf", width = 6, height = 6)
plotEffects(effects4plot = data4plot$alleffs, xmargs = c(0.1, 0.3), label_size = 1.5,
            sortlabs = 1:11, title_text = "")
dev.off()


