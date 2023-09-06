runMRnet <- function(dataExpression, writeCsv=FALSE, pathOut="."){
  mim <- minet::build.mim(dataExpression, estimator = "spearman")
  net <- minet::minet(mim ,method="mrnet", estimator="spearman")
  if(writeCsv){
    #fileName <- paste0(pathOut,"mrnet_predicted_original.csv")
    #writeNetworkCsv(netMrnet, fileName)
    fileR = paste0(pathOut,"mrnet_predicted_original.RData")
    writeRData(net, fileR)
  }
  return(net)
}