removeGeneGeneNode <- function(tfList, netOrig, writeData = FALSE, pathOut = "."){
  message(paste((sum(netOrig == 1)),"Node found before TF filter"))
  valid_names <- c(tfList[[1]])
  net <- as.matrix(netOrig)
  
  for (row_name in rownames(net)) {
    for (col_name in colnames(net)) {
      if (net[row_name, col_name] == 1 && !(row_name %in% valid_names) && !(col_name %in% valid_names)) {
        net[row_name, col_name] <- 0
      }
    }
  }
  message(paste((sum(net == 1)),"Node found after TF filter"))
  if(writeData){
    message("Writting data...")
    fileName <- paste0(pathOut,"filtred_geneGene.csv")
    writeNetworkCsv(net, fileName)
    fileR = paste0(pathOut,"filtred_geneGene.RData")
    writeRData(net, fileR)
    remove(fileR)
  }
  remove(fileName)
  return(net)
}

