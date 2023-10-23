runGENIE3 <- function(dataExpression, writeData=FALSE, pathOut="."){
  set.seed(123)
  net <-GENIE3(t(dataExpression), nCores = 4, verbose = TRUE)
  if(writeData){
    #fileName <- paste0(pathOut,"aracne_predicted_original.csv")
    #writeNetworkCsv(net, fileName)
    fileR = paste0(pathOut,"aracne_predicted_original.RData")
    writeRData(net, fileR)
  }
  return(net)
}