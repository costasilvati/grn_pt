runC3net <- function(dataExpression, writeCsv=FALSE, pathOut="."){
  net <- c3net::c3net(c3net::makemim(t(dataExpression)))
  if(writeCsv){
    #file <- paste0(pathOut, "c3net_predicted_original.csv")
    #writeNetworkCsv(net, file)
    fileR = paste0(pathOut,"c3Net_predicted_original.RData")
    writeRData(net, fileR)
  }
  return(net)
}