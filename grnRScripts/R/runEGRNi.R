runEGRNi <- function(dataExpression, tfList='', writeCsv=FALSE, pathOut="."){

  net <- CRN(dataExpression)
  if(writeCsv){
    fileR = paste0(pathOut,"EGRNi_predicted_original.RData")
    #save(net, file = fileR)
    writeRData(net, fileR)
  }
  return(net)
}
