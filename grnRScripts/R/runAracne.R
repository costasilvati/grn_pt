runAracne <- function(dataExpression, writeCsv=FALSE, pathOut="."){
  mim <- minet::build.mim(dataExpression)
  net <- minet::aracne(mim)
  if(writeCsv){
    #fileName <- paste0(pathOut,"aracne_predicted_original.csv")
    #writeNetworkCsv(net, fileName)
    fileR = paste0(pathOut,"aracne_predicted_original.RData")
    writeRData(net, fileR)
  }
  return(net)
}