runClr <- function(dataExpression,writeCsv = FALSE, pathOut = "."){
  netClr <- minet::minet(dataset = dataExpression, method = "clr")
  if(writeCsv){
    fileName <- paste0(pathOut,"clr_predicted_original.csv")
    writeNetworkCsv(netMinetClr, fileName)
  }
  return(netClr)
}