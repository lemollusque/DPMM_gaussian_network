## ---------------------------------------------------------------------------------------------------------------------------------------------------------
## load libraries
## ---------------------------------------------------------------------------------------------------------------------------------------------------------
require(BiDAG) ## for DAG sampling
require(Bestie) ## for intervention effects

require("excel.link")
require(readxl)
require(corrplot)
require(dplyr)
require(tidyverse)
require(magrittr)
require(foreach)
require(doParallel)
require(parallelly)
require(DiagrammeR)
require(DiagrammeRsvg) ## for exporting svg for plotting to file
require(rsvg) ## for converting svg to png

# load functions script
source("functions.R")

# load env variables
require(dotenv)
load_dot_env(file = ".env")

## ---------------------------------------------------------------------------------------------------------------------------------------------------------
## Load and prepare the data
## ---------------------------------------------------------------------------------------------------------------------------------------------------------
# import data
year = "2007"
complete_data = importData(year=year)

# weights
weights = complete_data$weight_score

inputData = preprocessData(data=complete_data)

# dataname
dataname <- paste0("data_", year)
## ---------------------------------------------------------------------------------------------------------------------------------------------------------
## DAG sampling
## ---------------------------------------------------------------------------------------------------------------------------------------------------------
nDAGs <- 100
batch <- 101:102
nSeeds <- length(batch)
labels4plot <- colnames(inputData)
## equivalent to names(ncs4analysis)
nNodes <- length(labels4plot)

## Set edges blacklist
blacklist = matrix(0, length(colnames(inputData)), length(colnames(inputData)))
## childhood sexual abuse and bullying only admit incoming edges from sex and from each other
blacklist[,c(which(colnames(inputData)=="Sex_abuse_before_16"), which(colnames(inputData)=="Bully"))] = 1
blacklist[which(colnames(inputData)=="Sex"), c(which(colnames(inputData)=="Sex_abuse_before_16"), which(colnames(inputData)=="Bully"))] = 0
blacklist[which(colnames(inputData)=="Bully"), which(colnames(inputData)=="Sex_abuse_before_16")] = 0
blacklist[which(colnames(inputData)=="Sex_abuse_before_16"), which(colnames(inputData)=="Bully")] = 0
## adult sexual abuse only admits incoming edges from sex, bullying and earlier abuse
blacklist[, which(colnames(inputData)=="Sex_abuse_since_16")] = 1
blacklist[which(colnames(inputData)=="Sex"), which(colnames(inputData)=="Sex_abuse_since_16")] = 0
blacklist[which(colnames(inputData)=="Bully"), which(colnames(inputData)=="Sex_abuse_since_16")] = 0
blacklist[which(colnames(inputData)=="Sex_abuse_before_16"), which(colnames(inputData)=="Sex_abuse_since_16")] = 0
## while sex only admits outgoing edges since symptoms cannot cause it.
blacklist[, which(colnames(inputData)=="Sex")] = 1




# use parallelization for sampling
numCores <- min(length(batch), parallelly::availableCores())
cl <- makeCluster(numCores)
registerDoParallel(cl)
# sampling loop
foreach(seednumber=batch) %dopar% {
  timing <- proc.time()
  print(paste("Seed is", seednumber))
  sampleDAGs(inData=inputData,
             scoretype="bde",
             bdepar = list(edgepf = 1, chi = 1),
             nDigraphs=nDAGs,
             seed=seednumber,
             dname=dataname,
             weightvector=weights,
             edgeBlacklist=blacklist,
             bgnodes=c(which(colnames(inputData)=="Sex")))
  computeEffects(n=nNodes, ## careful if files do not exist :)
                 seed=seednumber,
                 dname=dataname)
  print(proc.time() - timing)
}
stopCluster(cl)

## ---------------------------------------------------------------------------------------------------------------------------------------------------------
## Draw DAG
## ---------------------------------------------------------------------------------------------------------------------------------------------------------
data4plot <- loadsamples(seeds=batch, nn=nNodes, dname=dataname)
allDAGs <- data4plot$alldigraphs
allEffects <- data4plot$alleffs


# plot the DAG
graph2plot <- dagviz(allDAGs,
                     grouped_nodes = list(c(which(colnames(inputData)=="Sex_abuse_before_16"), which(colnames(inputData)=="Bully"))),
                     style_mat = matrix(1, nrow=length(labels4plot), ncol=length(labels4plot)))
displayDAG(g2plot = graph2plot, figname=paste0("DAG_", dataname, ".png"))

# prepare labels for effects plot
labels_in_order = c("Sex",
                    "Bully",
                    "Sex_abuse_before_16",
                    "Sex_abuse_since_16",
                    "Paranoia",
                    "Worry",
                    "Mood_instability",
                    "Fatigue",
                    "Phobias",
                    "Depression",
                    "Irritability",
                    "Concentration",
                    "Cannabis",
                    "Compulsions",
                    "Hallucinations",
                    "Obsessions",
                    "Panic",
                    "Physical_worry",
                    "Sleep",
                    "Somatic")

labels_order = c()
for (l in labels_in_order){
  labels_order = c(labels_order, which(colnames(inputData)==l))
}

# plot the effects
plotEffects(effects4plot = allEffects, # effects4plot$allarray,
            labs = labels4plot, 
            sortlabs = labels_order, 
            figname=paste0("causal_effects_", dataname, ".png"), 
            title_text = "Distributions of Downstream Causal Effects")

