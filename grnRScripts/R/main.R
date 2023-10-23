setwd("grnRScripts/R/")
source("loadAll.R")
loadAll()

#---------

pathFile = "C:\Users\jujuc\OneDrive\Doutorado\Redes_PT\Network3_expression_data.tsv"
pathFileTf = "C:\Users\jujuc\OneDrive\Doutorado\Redes_PT\Network3_transcription_factors.tsv"
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

networks <- inferenceNetwork(expressionData, TRUE, pathOut) # R tools
# load("/Volumes/SD128/GRN_PT/netDream5Net3/networkInference/listAllNetworks_predicted_predicted.RData")
# networks <- data
# remove(data)
#------ Pyhton Inference Tools --------
pathPythonTools <- "/Users/julianacostasilva/Google\ Drive/Meu\ Drive/RNA-Seq\ Herbas/GRNs/resultsDream5Net3/"
GENIE3 <- read.csv(paste0(pathPythonTools, "GENIE3_predicted_original.csv"), 
                   header = TRUE, 
                   row.names = "X")
networks[["GENIE3"]] <- as.matrix(GENIE3)
# Verificar execução do ENNET
EnNET <- read.csv(paste0(pathPythonTools, "ennet_predict_original.csv"), 
                  header = TRUE, 
                  row.names = "X")
networks[["EnNET"]] <- as.matrix(EnNET)
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

#------------- Read gold Standard network and TF List ----------
tfList <- read_csv(pathFileTf, col_names = FALSE)

goldStndardMatrix <- listToMatrix(pathListFile= pathFileGoldSt,
                                  colRowNames=c(colnames(networks[[1]])), 
                                  dimension = 4511, 
                                  headerExsit=FALSE,
                                  separator = "\t",
                                  writeRData = TRUE, 
                                  pathOut = paste0(pathOut,"goldStandard"))
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
write.csv(results, file = paste0(pathOut,"allNetworksMetricas_0910.csv"))
# GeCoNet 0 na remoção GeneGeneConnection
#minThreshold <- min(net)
#maxThreshold <- max(net)
#stepIncrement = (maxThreshold-minThreshold)/nSteps
writeRData(results, paste0(pathOut,"metricasMultipleThresholds_2709.RData"))
write.csv(results, file = paste0(pathOut,"metricasMultipleThresholds_2609.csv"))
remove(namesRow, num_linhas, num_colunas, pattern, files, pattern2, countLines)
#------ Consenso
toolNames <- names(networks)
networksBest <- networks
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
    paste0(pathOut,"multipleThresholds/geneGeneConnection/BC3NET_0.1_filtred_geneGene.RData"),
    paste0(pathOut,"multipleThresholds/geneGeneConnection/C3NET_0.1781_filtred_geneGene.RData"),
    paste0(pathOut,"multipleThresholds/geneGeneConnection/ARACNE_0.1_filtred_geneGene.RData"),
    paste0(pathOut,"multipleThresholds/geneGeneConnection/CLR_0.2_filtred_geneGene.RData"),
    paste0(pathOut,"multipleThresholds/geneGeneConnection/mrnet_0.1_filtred_geneGene.RData"),
    paste0(pathOut,"multipleThresholds/geneGeneConnection/mrnetB_0.3_filtred_geneGene.RData"),
    paste0(pathOut,"multipleThresholds/geneGeneConnection/GENIE3_0.0511_filtred_geneGene.RData"),
    paste0(pathOut,"multipleThresholds/geneGeneConnection/EnNET_0.0_filtred_geneGene.RData"),
    paste0(pathOut,"multipleThresholds/geneGeneConnection/GeCoNet_0_filtred_geneGene.RData")
)
i <- 1 #contador de toolNames
for(f in files){
  load(f)
  networksBest[[i]] <- data
  remove(data)
  i <- i + 1
}
remove(files, f, i)
# ------------- Gerar rede consenso ------------------
netConsenso <- consensus(networks = networksBest, goldStndardMatrix, namesCol, writeData = TRUE)

write.csv(netConsenso$net, file=paste0(pathOut,"consensusNet.csv"))

netConsenso$result$threshold[1] <- 0.1

#----- ROC Curve ------------------
library(PRROC)
toolNames <- names(networks)

for (tool in toolNames) {
  linesTool <- grepl(tool, row.names(results))
  dataTool <- results[linesTool, ]
  pr_result <- pr.curve(scores.class0 = as.matrix(goldStndardMatrix), scores.class1 = as.matrix(dataTool))
  print(paste(tool," AUPR ",pr_resu))
  roc <- pr.curve(dataTool$precision,dataTool$sensitivity, curve=TRUE)
  plot(roc, main=tool)
}

# remove(BC3NET, data, dataTool, edges_gold, edges_net, EnNET, GeCoNet, GENIE3, linhas_colunas_top_100, netPy, netT, resultsTemp, triangular_superior)
# remove(aupr_values, countLines, fn, fp, namesCol, top_100000_valores_mutua, tp, valores_mutua, valores_mutua_ordenados)

# AUPR
install.packages("pROC")
library(pROC)

# Calcular a curva ROC
roc_obj <- roc(as.vector(nets$ARACNE), as.vector(goldStndardMatrix))

# Calcular a AUPR
aupr <- auc(roc_obj, method = "PR")

# Exibir o valor da AUPR
print(paste("AUPR:", round(aupr, 3)))
#--- TESTES ----
# Carregue as bibliotecas necessárias
library(pROC)
library(readr)
library(dplyr)
gold <- read_csv("/Volumes/Macintosh\ HD/Users/julianacostasilva/Google\ Drive/Meu\ Drive/RNA-Seq\ Herbas/GRNs/resultsDream5Net3/bestThresholds/GENIE3_0.0510852_filtred_geneGene.csv")

gold <- as.data.frame(gold)

# Defina as colunas do dataframe como nomes de linha
rownames(gold) <- gold$...1

# Remova a primeira coluna se necessário
gold <- gold[, -1]
gold <- as.vector(unlist(gold))

# Converta seu dataframe 'df' para um vetor
df <- as.vector(unlist(df))

# Encontre os índices onde df é igual a 1
indices <- which(df == 1)

# Extraia os valores correspondentes aos índices
pred <- as.numeric(df[indices])
real <- as.numeric(gold[indices])

# Calcule a precisão e a revocação
precision <- cumsum(real) / (1:length(real))
recall <- cumsum(real) / sum(real)

roc_obj <- roc(real, pred)
aupr <- auc(roc_obj, method = "sens_spec")
aupr <- aup

