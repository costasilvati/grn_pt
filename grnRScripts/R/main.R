setwd("grnRScripts/R/")
source("loadAll.R")
loadAll()

#---------

pathFile = "C:\Users\jujuc\OneDrive\Doutorado\Redes_PT\Network3_expression_data.tsv"
pathFileTf = "C:\Users\jujuc\OneDrive\Doutorado\Redes_PT\Network3_transcription_factors.tsv"
pathFileGoldSt <- "/Volumes/SD128/GRN_PT/DREAM5_NetworkInference_GoldStandard/DREAM5_NetworkInference_GoldStandard_Network3.tsv"
pathOut = "/Volumes/SD128/GRN_PT/netDream5Net3/"
#----------- Lendo os dados de expressão e goldStandard -------------
expressionData <- read.delim(pathFile)
writeRData(expressionData, paste0(pathOut,"expressionData.RData"))
tfList <- read_csv(pathFileTf, col_names = FALSE)

# ------ Lendo lista de nós Gold Standard e convertendo em matrix ---------------------
goldStandard <- read.delim(pathFileGoldSt, 
                           header = FALSE) # edge list
writeRData(goldStandard, paste0(pathOut,"goldStandard.RData"))
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


