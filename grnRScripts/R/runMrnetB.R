runMrnetB <- function(dataExpression, write=FALSE, pathOut="."){
  mim <- minet::build.mim(dataExpression, estimator = "spearman")
  net <- minet::minet(mim ,method="mrnetb", estimator="spearman")
  if(write){
    #fileR = paste0(pathOut,"mrnetB_predicted_original.RData")
    #writeRData(net, fileR)
    fileR = paste0(pathOut,"mrNetB_predicted_original.RData")
    writeRData(net, fileR)
  }
  return(net)
}