runEnnet <- function(dataExpression,tf = NULL, writeCsv=FALSE, pathOut="."){
  if(is.null(tf)){
    netEnnet <- ennet::ennet(E=dataExpression)
  }else{
    netEnnet <- ennet::ennet(E=dataExpression, Tf=tf)
  }
  if(writeCsv){
    fileName <- paste0(pathOut,"ennet_predicted_original.csv")
    writeNetworkCsv(netEnnet, fileName)
  }
  return(netEnnet)
}