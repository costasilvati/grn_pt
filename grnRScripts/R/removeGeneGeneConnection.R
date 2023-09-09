# Pattern output name: pathOut + _ + tool + threshold + filtred_geneGene.RData
removeGeneGeneConnection <- function(tfList, netOrig, writeData = FALSE, pathOut = "."){
  nos_na_lista <- colnames(netOrig) %in% tfList
  netOrig[!nos_na_lista, !nos_na_lista] <- 0
  message(paste((sum(netOrig == 1))," Connection found after TF filter"))
  if(writeData){
    message("Writting data...")
    # fileName <- paste0(pathOut,"filtred_geneGene.csv")
    # writeNetworkCsv(netOrig, fileName)
    fileR = paste0(pathOut,"filtred_geneGene.RData")
    writeRData(netOrig, fileR)
    remove(fileR)
  }
  #return(net)
}

removeGeneGeneNode2 <- function(tfList, netOrig, writeData = FALSE, pathOut = ".") {
  valid_names <- c(tfList[[1]])
  net <- as.matrix(netOrig)
  
  net <- sapply(rownames(net), function(row_name) {
    sapply(colnames(net), function(col_name) {
      if (net[row_name, col_name] == 1 && !(row_name %in% valid_names) && !(col_name %in% valid_names)) {
        return(0)
      } else {
        return(net[row_name, col_name])
      }
    })
  })
  
  message(paste((sum(net == 1)), "Node found after TF filter"))
  
  if (writeData) {
    message("Writting data...")
    # fileName <- paste0(pathOut, "filtred_geneGene.csv")
    # writeNetworkCsv(net, fileName)
    fileR <- paste0(pathOut, "filtred_geneGene.RData")
    writeRData(net, fileR)
    remove(fileR)
  }
  
  return(net)
}

