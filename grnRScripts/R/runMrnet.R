runMRnet <- function(dataExpression, writeCsv=FALSE, pathOut="."){
  netMRnet <- minet::minet(dataset = dataExpression, method = "mrnet")
  if(writeCsv){
    fileName <- paste0(pathOut,"mrnet_predicted_original.csv")
    writeNetworkCsv(netMRnet, fileName)
  }
  return(netMRnet)
}