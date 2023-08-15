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
#networks <- inferenceNetwork(expressionData, TRUE, pathOut)
#--------  ou importar RData
load("/Volumes/SD128/GRN_PT/netDream5Net3/listAllNetworks_predicted_predicted.RData")
networks <- data
remove(data)

names_net <- names(networks)
cont <- 1

networks_filtred <- list(names(networks))

for (net in networks) {
  message(names_net[cont])
  thresholdNetwork <- quantileByData(net)
  networks_filtred[[cont]] <- edgeMatrixByTreshold(network=net, 
                                      threshold= thresholdNetwork, 
                                      writeData = TRUE, 
                                      pathOut = paste(pathOut, names_net[cont])
                                      )
  cont <- cont + 1
}
names(networks_filtred) <- names_net 

networks_filtred_tf <- list(names_net)

cont2 <- 1
for (netF in networks_filtred) {
  message(names_net[cont2])
  net_only_TfGene <- removeGeneGeneNode(tfList,
                                        netF,
                                        TRUE, 
                                        paste0(pathOut,names_net[cont2]))
  networks_filtred_tf[[cont2]] <- net_only_TfGene
  cont2 <- cont2 + 1
}
names(networks_filtred_tf) <- names_net 
remove(cont)
remove(cont2)
remove(names_net
       )
remove(net)
remove(netF)
remove(networks_filtred)
remove(net_only_TfGene)


#---------------------- GENIE3--------
#netGenie3 <- read.csv("/Users/julianacostasilva/Google\ Drive/.shortcut-targets-by-id/1aap66FUWBqx6ORLHbTivskehalbbaq8R/RNA-Seq\ Herbas/GRNs/resultsDream5Net3/GENIE3_predicted_original.csv")
#netGenie3_treshold <- quantileByData(netGenie3)
#netGenie3_edgeList <- edgeListByTreshold(netGenie3, treshold = 0.01)
#netGenie3_edgeListTF <- removeGeneGeneNode(tfList, 
#                                           netGenie3_edgeList,
#                                           TRUE,
#                                           paste0(pathOut,"netGenie3_"))

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
row.names(goldSatndardMatrix) <- row.names(netBc3)
colnames(goldSatndardMatrix) <- row.names(netBc3)
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
#------------ Compare net---------



