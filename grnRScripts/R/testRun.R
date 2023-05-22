setwd("grnRScripts/R/")

source("readExpressionData.R")
source("runAracne.R")
source("runBc3Net.R")
source("runC3net.R")
source("runClr.R")
source("runEnnet.R")
source("runMrnet.R")
source("runMrnetB.R")
source("writeNetworkCsv.R")
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
dream5Net3 <- read.delim(pathFile)
#gera matriz transposta
dream5Net3T <- t(dream5Net3)
dream5Net3Tf <- read.delim(pathFileTf)
#executa bc3net
netBc3 <- runBc3Net(dream5Net3T,TRUE,pathOut)
#executa C3net
netC3 <- runC3net(dream5Net3T, TRUE, pathOut)
#executa CLR - não transposto
netCrl <- runClr(dream5Net3, TRUE, pathOut)
# executa mrnet - não transposto
netMrnet <- runMrnet(dream5Net3, TRUE, pathOut)
# executa MRNETB
netMrnetB <- runMrnet(dream5Net3T, TRUE, pathOut)
# executa ENNET
ennetNet <- runEnnet(dream5Net3T,dream5Net3Tf, TRUE, pathOut)

# read list edges gold standard
goldNet3 <- read.delim(pathFileGoldSt, header = FALSE)
adjGold <- matrix(data = 0, nrow = 4511, ncol = 4511)
row.names(adjGold) <- row.names(mim)
colnames(adjGold) <- row.names(mim)
node1 <- goldNet3[1]
node2 <- goldNet3[2]
edge <- goldNet3[3]
line <- 1
for (n1 in node1[,1]) {
  adjGold[n1,node2[line, 1]] <- edge[line, 1]
  line <- line +1
}


