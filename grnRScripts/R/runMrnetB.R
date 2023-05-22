runMrnetB <- function(dataExpression, writeCsv=FALSE, pathOut="."){
  mim <- minet::build.mim(dataExpression,estimator="spearman")
  netMRNetB <- minet::mrnetb(mim)
  if(writeCsv){
    fileName <- paste0(pathOut,"mrnetB_predicted_original.csv")
    writeNetworkCsv(netMRNetB, fileName)
  }
  return(netMRNetB)
}