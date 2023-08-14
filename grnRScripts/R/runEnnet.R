runEnnet <- function(dataExpression,tf = NULL, writeCsv=FALSE, pathOut="."){
  if(is.null(tf)){
    net <- ennet::ennet(E=dataExpression)
  }else{
    net <- ennet::ennet(E=dataExpression, Tf=tf)
  }
  if(writeCsv){
    #fileName <- paste0(pathOut,"ennet_predicted_original.csv")
    #writeNetworkCsv(netEnnet, fileName)
    fileR = paste0(pathOut,"ennet_predicted_original.RData")
    writeRData(netMRNetB, fileR)
  }
  return(net)
}