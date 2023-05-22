runAracne <- function(dataExpression, writeCsv=FALSE, pathOut="."){
  mim <- minet::build.mim(dataExpression,estimator="spearman")
  netAracne <- minet::aracne(mim)
  if(writeCsv){
    fileName <- paste0(pathOut,"aracne_predicted_original.csv")
    writeNetworkCsv(netAracne, fileName)
  }
  return(netAracne)
}