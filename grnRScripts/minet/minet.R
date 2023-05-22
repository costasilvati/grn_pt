source("readExpressionData.R")
library(minet)
BiocManager::install("Rgraphviz")
library(Rgraphviz)

# ARACNE
pathFile <- "/Volumes/SD128/GRN_PT/DREAM5_NetworkInferenceChallenge_Data/Network3_Data/Network3_expression_data.tsv"
dataExp <- read.delim(pathFile)

mim <- minet::build.mim(data,estimator="spearman")
netAracne <- minet::aracne(mim)
#plot(as(netAracne ,"graphNEL"))

#MRNetB
netMRNetB <- minet::mrnetb(mim)

#minet method aracne
netMinetAracne <- minet::minet(dataset = data, method = "aracne")

#mrnet
netMinetMRnet <- minet::minet(dataset = data, method = "mrnet")

#clr
netMinetClr <- minet::minet(dataset = data, method = "clr")

# Converter tsv GoldStandard em matriz