runC3net <- function(dataExpression, writeCsv=FALSE, pathOut="."){
  net <- c3net::c3net(dataExpression)
  if(writeCsv){
    file <- paste0(pathOut, "c3net_predicted_original.csv")
    writeNetworkCsv(net, file)
  }
  return(net)
}