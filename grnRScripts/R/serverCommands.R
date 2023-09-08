# setwd("grnRScripts/R/")
library(dplyr)
source("loadAll.R")
loadAll()
load("../data/listAllNetworks_predicted_predicted.RData")
networks <- data
remove(data)
load("../data/goldStandard.RData")
goldStandard <- data
remove(data)
load("../data/goldStandardTranscFactor.RData")
tfList <- data
remove(data)

pathOut <- "../../../../grn/"
expressionData <- read.delim("../../../../grn/Network3_expression_data.tsv")

networks <- inferenceNetwork(expressionData, TRUE, pathOut)

# Gerar vários filtros
names_net <- names(networks)
netThresolds <- list()
netThresolds <- names(names_net)
i <- 1
for (net in networks) {
  minThreshold <- min(net)
  nameList <- names_net[i]
  namesT <- NULL
  t <- 1
  tempList <- list()
  while (minThreshold <= 1.0) {
    message(paste(nameList," filter by threshold: ",minThreshold))
    netF <- edgeMatrixByTreshold(net,minThreshold,writeData = TRUE)
    message(paste(nameList," removing gene-gene connections"))
    tempList[paste0(nameList,minThreshold)] <- removeGeneGeneNode2(tfList,netF,FALSE)
    minThreshold <- minThreshold + 0.1
    t <- t+1
  }
  netThresolds[[nameList]] <- tempList
  i <- i+1
}
remove(i,t,tempList, minThreshold, namesT, netF, net)