setwd("grnRScripts/R/")

source("readExpressionData.R")
source("runAracne.R")
source("runBc3Net.R")
source("runC3net.R")
source("runClr.R")
source("runEnnet.R")
source("runMrnet.R")
source("runMrnetB.R")
source("writeRData.R")
source("writeNetworkCsv.R")
source("quantileByData.R")
source("filterNetworkByTreshold.R")
source("removeGeneGeneNode.R")
library('bc3net')
library('c3net')
library('wrMisc')
library(readr)
library(dplyr)
library(minet)
library(ennet)
library(netdiffuseR)
library(tidyverse)
library(tibble)

#---------

pathFile = "/Volumes/SD128/GRN_PT/DREAM5_NetworkInferenceChallenge_Data/Network3_Data/Network3_expression_data.tsv"
pathFileTf = "/Volumes/SD128/GRN_PT/DREAM5_NetworkInferenceChallenge_Data/Network3_Data/Network3_transcription_factors.tsv"
pathFileGoldSt <- "/Volumes/SD128/GRN_PT/DREAM5_NetworkInference_GoldStandard/DREAM5_NetworkInference_GoldStandard_Network3.tsv"
pathOut = "/Volumes/SD128/GRN_PT/netDream5Net3/"
#----------- Lendo os dados de expressão e goldStandard --------------------------
expressionData <- read.delim(pathFile)
# writeRData(expressionData, paste0(pathOut,"expressionData.RData")) #
# -- gera matriz transposta
expressionDataTransposto <- t(expressionData)
# writeRData(expressionData, paste0(pathOut,"expressionDataTransposto.RData"))
#----------- Inferindo e filtando redes --------------------------
# ------ bc3net - MI qu-------
netBc3 <- runBc3Net(expressionDataTransposto,
                    TRUE,
                    pathOut)
# ------ bc3net - Calcula e filtra por quartil 4 ------------------------
netBc3treshold = quantileByData(netBc3)
netBc3_edgeListQuantile <- edgeListByTreshold(netBc3, treshold = netBc3treshold)
netBc3_edgeList <- edgeListByTreshold(netBc3, treshold = 0.01)
# --- bc3net - Remove conexão gene - gene  --------------
netBc3_edgeListTF <- removeGeneGeneNode(tfList, 
                                        netBc3_edgeListQauntile,
                                        TRUE,
                                        paste0(pathOut,"netBc3_"))
# ------ C3net - MI -------
netC3 <- runC3net(expressionDataTransposto, 
                  TRUE, 
                  pathOut)
# ------ C3net - Calcula e filtra por quartil 4 ------------------------
netC3Trshold <- quantileByData(netC3)
netC3_edgeListQauntile <- edgeListByTreshold(netC3, treshold = netC3Trshold)
netC3_edgeList <- edgeListByTreshold(netC3, treshold = 0.01)

# --- C3net - Remove conexão gene - gene  --------------
netC3_edgeListTF <- removeGeneGeneNode(tfList, 
                                        netC3_edgeList,
                                        TRUE,
                                        paste0(pathOut,"netC3_"))
# -- executa Aracne - transposto?
netAracne <- runAracne(expressionData, TRUE, pathOut)
netAracneTreshold <- quantileByData(netAracne)
netAracne_edgeListQauntile <- edgeListByTreshold(netAracne, treshold = netAracneTreshold)
netAracne_edgeList <- edgeListByTreshold(netAracne, treshold = 0.01)

# -- executa CLR - não transposto
netClr <- runClr(expressionData,  TRUE, pathOut)
netClrTreshold <- quantileByData(netClr)
netClr_edgeListQuantile <- edgeListByTreshold(netClr, treshold = netClrTreshold)
netClr_edgeList <- edgeListByTreshold(netClr, treshold = 0.01)
netClr_edgeListTF <- removeGeneGeneNode(tfList, 
                                           netClr_edgeList,
                                           TRUE,
                                           paste0(pathOut,"netClr_"))
#---------------------- GENIE3--------
netGenie3 <- read.csv("/Users/julianacostasilva/Google\ Drive/.shortcut-targets-by-id/1aap66FUWBqx6ORLHbTivskehalbbaq8R/RNA-Seq\ Herbas/GRNs/resultsDream5Net3/GENIE3_predicted_original.csv")
netGenie3_treshold <- quantileByData(netGenie3)
netGenie3_edgeList <- edgeListByTreshold(netGenie3, treshold = 0.01)
netGenie3_edgeListTF <- removeGeneGeneNode(tfList, 
                                           netGenie3_edgeList,
                                           TRUE,
                                           paste0(pathOut,"netGenie3_"))

# -- executa mrnet - não transposto
netMrnet <- runMRnet(expressionData, 
                     TRUE, 
                     pathOut)
netMrnetTreshold <- quantileByData(netMrnet)
netMrnet_edgeListQuantile <- edgeListByTreshold(netMrnet, treshold = netMrnetTreshold)
netMrnet_edgeList <- edgeListByTreshold(netMrnet, treshold = 0.01)
netMrnet_edgeListTF <- removeGeneGeneNode(tfList, 
                                        netMrnet_edgeList,
                                        TRUE,
                                        paste0(pathOut,"netMrnet_"))
# -- executa MRNETB
netMrnetB <- runMrnetB(expressionDataTransposto, TRUE, pathOut)
netMrnetBTreshold <- quantileByData(netMrnetB)
netMrnetB_edgeListQuantile <- edgeListByTreshold(netMrnetB,netMrnetBTreshold)
netMrnetB_edgeList <- edgeListByTreshold(netMrnetB,0.01)
netMrnetB_edgeListTF <- removeGeneGeneNode(tfList, 
                                          netMrnetB_edgeList,
                                          TRUE,
                                          paste0(pathOut,"netMrnetB_"))
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
node1 <- goldNet3[1]
node2 <- goldNet3[2]
edge <- goldNet3[3]
#  -- create edge in matrix
line <- 1
for (n1 in node1[,1]) {
  goldSatndardMatrix[n1,node2[line, 1]] <- edge[line, 1]
  line <- line +1
}
#------------ Compare net---------



