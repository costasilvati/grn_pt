library(dplyr)
removeGeneGeneNode <- function(transcFactorList, edgeList, writeData = FALSE, pathOut = "."){
  lines_out = c()
  out = 1
  for (l in 1:nrow(edgeList)) {
    if(edgeList[l,1] %in% transcFactorList | edgeList[l,2] %in% transcFactorList){
      message(paste("TF - Gene in line ", l))
    }else{
      lines_out[out] = l
      out = out+1
      message(paste("Gene - Gene in line ", l))
    }
  }
  edgeList <- edgeList %>% filter(!dplyr::row_number() %in% lines_out)
  if(writeData){
    fileName <- paste0(pathOut,"filtred_geneGene.csv")
    writeNetworkCsv(edgeList, fileName)
    fileR = paste0(pathOut,"filtred_geneGene.RData")
    writeRData(edgeList, fileR)
  }
  return(edgeList)
}