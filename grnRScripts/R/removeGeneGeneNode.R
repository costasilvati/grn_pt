removeGeneGeneNode <- function(transcFactorList, netMatrix, writeData = FALSE, pathOut = "."){
  lineNames <- row.names(netMatrix)
  columnNames <- row.names(netMatrix)
  tfList <- c(transcFactorList[[1]])
  nodeFound <- 0
  nodeRemoved <- 0
  for (line in lineNames) {
    for (column in columnNames) {
      if(netMatrix[line, column] == 1){ # se nó existe
        nodeFound <- nodeFound + 1
        if(!((line %in% tfList) | (column %in% tfList))){
          netMatrix[line, column] <- 0
          nodeRemoved <- nodeRemoved + 1
        }
      }
    }
  }
  msg <- paste("Summary: \n Node founds: ", nodeFound)
  message(paste(paste(msg, "\n Node removed: "),nodeRemoved))
  report <- list(nodeFound, nodeRemoved, netMatrix)
  if(writeData){
    fileName <- paste0(pathOut,"filtred_geneGene.csv")
    writeNetworkCsv(report, fileName)
    fileR = paste0(pathOut,"filtred_geneGene.RData")
    writeRData(report, fileR)
  }
  return(report)
}