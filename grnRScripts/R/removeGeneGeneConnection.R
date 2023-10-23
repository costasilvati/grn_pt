# Pattern output name: pathOut + _ + tool + threshold + filtred_geneGene.RData
removeGeneGeneConnection <- function(tfList, netOrig, writeData = FALSE, pathOut = "."){
  tfList <- c(tfList)
  nos_na_lista <- colnames(netOrig) %in% tfList$X1
  netOrig[!nos_na_lista, !nos_na_lista] <- 0
  cat(sum(netOrig == 1)," Connection found after TF filter\n")
  if(writeData){
    message("Writting data...")
    # fileName <- paste0(pathOut,"_filtred_geneGene.csv")
    # writeNetworkCsv(netOrig, fileName)
    fileR = paste0(pathOut,"_filtred_geneGene.RData")
    writeRData(netOrig, fileR)
    remove(fileR)
  }
  return(netOrig)
}

