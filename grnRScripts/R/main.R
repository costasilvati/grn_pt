setwd("grnRScripts/R/")

source("readExpressionData.R")
source("runAracne.R")
source("runBc3Net.R")
source("runC3net.R")
source("runClr.R")
source("runEnnet.R")
source("runMrnet.R")
source("runMrnetB.R")
source("inferenceNetwork.R")
source("writeRData.R")
source("writeNetworkCsv.R")
source("inferenceNetwork.R")
source("quantileByData.R")
source("removeGeneGeneNode.R")
source("edgeMatrixByTreshold.R")
library(bc3net)
library(c3net)
library(wrMisc)
library(Matrix)
library(readr)
library(dplyr)
library(minet)
library(ennet)
library(netdiffuseR)
library(tidyverse)
library(tibble)
library(minet)
library(igraph)
library(stringr)
library(corto)

#---------

pathFile = "/Volumes/SD128/GRN_PT/DREAM5_NetworkInferenceChallenge_Data/Network3_Data/Network3_expression_data.tsv"
pathFileTf = "/Volumes/SD128/GRN_PT/DREAM5_NetworkInferenceChallenge_Data/Network3_Data/Network3_transcription_factors.tsv"
pathFileGoldSt <- "/Volumes/SD128/GRN_PT/DREAM5_NetworkInference_GoldStandard/DREAM5_NetworkInference_GoldStandard_Network3.tsv"
pathOut = "/Volumes/SD128/GRN_PT/netDream5Net3/"
#----------- Lendo os dados de expressão e goldStandard --------------------------
expressionData <- read.delim(pathFile)
writeRData(expressionData, paste0(pathOut,"expressionData.RData")) #
tfList <- read_csv(pathFileTf, col_names = FALSE)
#-------- Inferir redes 
networks <- inferenceNetwork(expressionData, TRUE, pathOut)
#--------  ou importar RData
# load("/Volumes/SD128/GRN_PT/netDream5Net3/listAllNetworks_predicted_predicted.RData")
# networks <- data
# remove(data)

# Filtrar pelo quartil 75%
#networks_filtred <- filterListNetwork(networks, TRUE, pathOut)

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
remove(i,t,tempList, minThreshold, namesT)

#------- 
goldStandardTranscFactor <- read.delim(pathFileTf)
writeRData(goldStandardTranscFactor, paste0(pathOut,
                                  "goldStandardTranscFactor.RData"))

# ------ Lendo lista de nós Gold Standard e convertendo em matrix ---------------------
goldStandard <- read.delim(pathFileGoldSt, 
                           header = FALSE) # edge list

# -- create matriz to gold standard data
goldSatndardMatrix <- matrix(data = 0, 
                             nrow = 4511, 
                             ncol = 4511) # matrix
row.names(goldSatndardMatrix) <- row.names(networks$netBc3)
colnames(goldSatndardMatrix) <- row.names(networks$netBc3)
# -- get node list
node1 <- goldStandard[1]
node2 <- goldStandard[2]
edge <- goldStandard[3]
#  -- create edge in matrix
line <- 1
for (n1 in node1[,1]) {
  goldSatndardMatrix[n1,node2[line, 1]] <- edge[line, 1]
  line <- line +1
}

sum(goldSatndardMatrix == 1)
sum(goldStandard == 1)
remove(line, node1, node2, edge)
#------------ Compare net---------
confusion_matrix <- list(names_net)
contConfusion <- 1
for (net in networks_filtred) {
  confusion_matrix[[contConfusion]] <- compareNetworks(goldSatndardMatrix, net)
  contConfusion <- contConfusion + 1
}
names(confusion_matrix) <- names_net
fileR <- paste0(pathOut,"confusion_matrix.RData")
writeRData(confusion_matrix, fileR)
remove(fileR, contConfusion)


