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
library(readr)
library(minet)
library(ennet)
library(netdiffuseR)

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
# ------ bc3net -------
netBc3 <- runBc3Net(expressionDataTransposto,
                    TRUE,
                    pathOut)
# ------ bc3net - Calcula e filtra por quartil 4 ------------------------
netBc3_filtred <- filterNetworkByTreshold(netBc3, 
                                          TRUE, 
                                          TRUE, 
                                          "Bc3")
# --- bc3net - Remove conexão gene - gene  --------------
netBc3_filtred2 <- removeGeneGeneNode(goldStandardTranscFactor, 
                                          netBc3_filtred, 
                                          TRUE,
                                          paste0(pathOut,"netBc3_"))
# ------ C3net -------
netC3 <- runC3net(expressionDataTransposto, 
                  TRUE, 
                  pathOut)
# ------ C3net - Calcula e filtra por quartil 4 ------------------------
netC3_filtred <- filterNetworkByTreshold(netC3, 
                                         TRUE, 
                                         TRUE, 
                                         pathOut, 
                                         "C3")
# --- C3net - Remove conexão gene - gene  --------------
netC3_filtred_2 <- removeGeneGeneNode(goldStandardTranscFactor, 
                                      netC3_filtred,
                                      TRUE,
                                      paste0(pathOut,"netBc3_"))
# -- executa CLR - não transposto
netClr <- runClr(expressionData, TRUE, pathOut)
netClr_filtred <- filterNetworkByTreshold(netCrl, TRUE, TRUE, pathOut, "Clr")
# -- executa mrnet - não transposto
netMrnet <- runMRnet(expressionData, TRUE, pathOut)
netMrnet_filtred <- filterNetworkByTreshold(netMrnet, TRUE, TRUE, pathOut, "Mrnet")
# -- executa MRNETB
netMrnetB <- runMrnetB(expressionData, TRUE, pathOut)
netMrnetB_filtred <- filterNetworkByTreshold(netMrnetB, TRUE, TRUE, pathOut, "MrnetB")
# -- executa ENNET
# ennetNet <- runEnnet(dream5Net3T,dream5Net3Tf, TRUE, pathOut)

goldStandardTranscFactor <- read.delim(pathFileTf)
writeRData(expressionData, paste0(pathOut,
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

# Falta remover gene-gene


#--------  Executa comparação com goldStandard ----------------
netBc3_valid_gold <- minet::validate(netBc3_filtred, goldSatndardMatrix)
writeRData(netBc3_valid_gold, paste0(pathOut,"Bc3_valid_gold.RData"))

netC3_valid_gold <- minet::validate(netC3_filtred, goldSatndardMatrix)
writeRData(netC3_valid_gold, paste0(pathOut,"C3_valid_gold.RData"))

netClr_valid_gold <- minet::validate(netClr_filtred, goldSatndardMatrix)
writeRData(netClr_valid_gold, paste0(pathOut, "Clr_valid_gold.RData"))

netMrnet_valid_gold <- minet::validate(netMrnet_filtred, goldSatndardMatrix)
writeRData(netMrnet_valid_gold, paste0(pathOut, "Mrnet_valid_gold.RData"))

netMrnetB_valid_gold <- minet::validate(netMrnetB_filtred, goldSatndardMatrix)
writeRData(netMrnetB_valid_gold, paste0(pathOut, "MrnetB_valid_gold.RData"))

#--- PLOT
max(fscores(netBc3_valid_gold))
dev <- show.pr(netBc3_valid_gold, col="green", type="b")
dev <- show.pr(netC3_valid_gold, device=dev, col="blue", type="b")
show.pr(netClr_valid_gold, device=dev, col="red",type="b")
auc.pr(netClr_valid_gold)

