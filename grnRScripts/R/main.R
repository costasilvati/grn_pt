setwd("grnRScripts/R/")
source("loadAll.R")
loadAll()

#---------

pathFile = "/Volumes/SD128/GRN_PT/DREAM5_NetworkInferenceChallenge_Data/Network3_Data/Network3_expression_data.tsv"
  #"C:\Users\jujuc\OneDrive\Doutorado\Redes_PT\Network3_expression_data.tsv"
pathFileTf = "/Volumes/SD128/GRN_PT/DREAM5_NetworkInferenceChallenge_Data/Network3_Data/Network3_transcription_factors.tsv"
  # "C:\Users\jujuc\OneDrive\Doutorado\Redes_PT\Network3_transcription_factors.tsv"
pathFileGoldSt <- "/Volumes/SD128/GRN_PT/DREAM5_NetworkInference_GoldStandard/DREAM5_NetworkInference_GoldStandard_Network3.tsv"
pathOut = "/Volumes/SD128/GRN_PT/netDream5Net3/"

#----------- Read Expression Data -------------
expressionData <- read.delim(pathFile)
writeRData(expressionData, paste0(pathOut,"expressionData.RData"))
# Falta 
# CORTO: https://cran.r-project.org/web/packages/corto/index.html
# iRafNet: https://cran.r-project.org/web/packages/iRafNet/index.html
# fssemR: https://cran.r-project.org/web/packages/fssemR/index.html
# EGRNi: https://cran.r-project.org/web/packages/EGRNi/index.html
#------------- Read gold Standard network and TF List ----------
tfList <- read_csv(pathFileTf, col_names = FALSE)

goldStndardMatrix <- listToMatrix(pathListFile= pathFileGoldSt,
                                  colRowNames=c(colnames(networks[[1]])), 
                                  dimension = 4511, 
                                  headerExsit=FALSE,
                                  separator = "\t",
                                  writeRData = TRUE, 
                                  pathOut = paste0(pathOut,"goldStandard"))
#-------------- Inference --------------
networks <- inferenceNetwork(expressionData, tfList, TRUE, pathOut) # R tools
# load("/Volumes/SD128/GRN_PT/netDream5Net3/networks.RData")
# networks <- data
# remove(data)
#------ Pyhton Inference Tools --------
# pathPythonTools <- "/Users/julianacostasilva/Google\ Drive/Meu\ Drive/RNA-Seq\ Herbas/GRNs/resultsDream5Net3/"
GeCoNet <- listToMatrix(pathListFile = paste0(pathPythonTools, "GeCoNet_predicted.csv"),
             colRowNames=c(colnames(networks[[2]])), 
             dimension = 4511, 
             headerExsit=TRUE,
             separator = ",",
             writeRData = TRUE, 
             pathOut = paste0(pathOut,"GeCoNet"))

networks[["GeCoNet"]] <- GeCoNet
remove(GeCoNet, GENIE3, EnNET)
#load(paste0(pathOut,"GeCoNet_Matrix.RData"))
#networks[["GeCoNet"]] <- data
#remove(data)


# ----- write RData with pattern name: tool_threshold_filtred_geneGene.RData
# Crie um data.frame de tamanho definido preenchido com zeros
nSteps <- 10 # How different thresholds need generate?
namesCol <-c("tool","threshold","fp","fn","tp","tn","accuracy","recall","precision","specificity","f_score","fdr")
results <- data.frame(matrix(0, nrow =(length(networks)* nSteps), ncol = length(namesCol)))
colnames(results) <- namesCol

results <- createMultipleThresholds(networks = networks, 
                                    goldNet = as.matrix(goldStndardMatrix),
                                    tfList = tfList, 
                                    nSteps = nSteps, # (max-min/nSteps)
                                    namesCol = namesCol,
                                    writeData = TRUE,
                                    pathOut=pathOut)
write.csv(results, file = paste0(pathOut,"allNetworksMetricas_071123.csv"))
# GeCoNet 0 na remoção GeneGeneConnection

#------ Consenso
toolNames <- names(networks)
networksBest <- results[[2]]
results[is.na(results)] <- 0.0
for (tool in toolNames) {
  linesTool <- grepl(tool, results$tool)
  dataTool <- results[linesTool, ]
  best_fscore <- max(dataTool$f_score)
  indices_linhas <- which(results$f_score == best_fscore)
  cat(tool, "\t best f_score: ", best_fscore," \t threshold:",results$threshold[indices_linhas[1]]," \n")
}
remove(toolNames, tool, linesTool, dataTool, best_fscore, indices_linhas)
# Com base na lista impressa resultante das linhas 85-95 , carregue os arquivos das redes com o 
#  melhor f-score, o nome do arquivo é formado por [nomeFerramenta_]+[threshold]+[_filtred_geneGene.RData].
#  Exemplo: ARACNE_0.1_filtred_geneGene.RData
# A ordem dos arquivos carregados deve seguir a ordem dos names da lista networks
files<- c(
    paste0(pathOut,"multipleThresholds/geneGeneConnection/BC3NET_0_filtred_geneGene.RData"),
    paste0(pathOut,"multipleThresholds/geneGeneConnection/C3NET_0_filtred_geneGene.RData"),
    paste0(pathOut,"multipleThresholds/geneGeneConnection/ARACNE_0.1_filtred_geneGene.RData"),
    paste0(pathOut,"multipleThresholds/geneGeneConnection/CLR_0.1_filtred_geneGene.RData"),
    paste0(pathOut,"multipleThresholds/geneGeneConnection/mrnet_0.1_filtred_geneGene.RData"),
    paste0(pathOut,"multipleThresholds/geneGeneConnection/mrnetB_0.2_filtred_geneGene.RData"),
    paste0(pathOut,"multipleThresholds/geneGeneConnection/GENIE3_0.0511_filtred_geneGene.RData"),
    paste0(pathOut,"multipleThresholds/geneGeneConnection/EnNET_0_filtred_geneGene.RData"),
    paste0(pathOut,"multipleThresholds/geneGeneConnection/GeCoNet_0_filtred_geneGene.RData"),
    paste0(pathOut,"multipleThresholds/geneGeneConnection/Corto_0.4546_filtred_geneGene.RData")
)
i <- 1 #contador de toolNames
networksBest <- networks
namesNetworks <- names(networksBest)
for(f in files){
  cat(namesNetworks[i])
  load(f)
  networksBest[[i]] <- data
  writeNetworkCsv(networksBest[[i]], paste0(pathOut,namesNetworks[i],"_bestThreshold.csv"))
  remove(data)
  i <- i + 1
}
remove(files, f, i)
# ------------- Gerar rede consenso ------------------
namesCol <-c("number of tool","threshold","fp","fn","tp","tn","accuracy","recall","precision","specificity","f_score","fdr")
resultsconsenso <- data.frame(matrix(0, nrow =length(networksBest), ncol = length(namesCol)))
colnames(resultsconsenso) <- namesCol

for (i in 1:length(networksBest)) {
  netConsenso <- consensus(networks = networksBest, 
                              gold =  goldStndardMatrix, 
                              namesCol = namesCol,
                              minConsenso = i, 
                              writeData = FALSE)
  resultsconsenso[i,] <- netConsenso$result
}


write.csv(netConsenso$net, file=paste0(pathOut,"consensusNet.csv"))

#------ Plot Venn e upset ------ 