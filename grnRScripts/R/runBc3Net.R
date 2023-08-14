runBc3Net <- function(dataExpression, writeCsv=FALSE, pathOut="."){
  net <-bc3net::bc3net(dataExpression, igraph=FALSE)
  if(writeCsv){
    #file = paste0(pathOut,"bc3net_predicted_original.csv")
    #writeNetworkCsv(net, file)
    fileR = paste0(pathOut,"bc3Net_predicted_original.RData")
    writeRData(net, fileR)
  }
  return(net)
}