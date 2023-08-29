inferenceNetwork <- function(expressionData, writeData = FALSE, pathOut = "."){
  # -- gera matriz transposta
  expressionDataTransposto <- t(expressionData)
  colnames(expressionDataTransposto) <- rep(1:length(row.names(expressionData)))
  message("Executing NetBc3...")
  netBc3 <- runBc3Net(expressionDataTransposto,writeData, pathOut)
  message("Executing NetC3...")
  netC3 <- runC3net(expressionDataTransposto, writeData, pathOut)
  message("Executing ARACNE...")
  netAracne <- runAracne(expressionDataTransposto, writeData, pathOut)
  message("Executing CLR...")
  netClr <- runClr(expressionDataTransposto, writeData, pathOut)
  message("Executing MRNet...")
  netMrnet <- runMRnet(expressionDataTransposto, writeData, pathOut)
  message("Executing MRNetB...")
  netMrnetB <- runMrnetB(expressionDataTransposto, writeData, pathOut)
  networks <- list(netBc3, netC3, netAracne, netClr, netMrnet, netMrnetB)
  names(networks) <- c("netBc3", "netC3", "netAracne", "netClr", "netMrnet", "netMrnetB")
  
  if(writeData){
    message("Writting data...")
    # fileName <- paste0(pathOut,"listAllNetworks_predicted.csv")
    # writeNetworkCsv(networks, fileName)
    fileR = paste0(pathOut,"listAllNetworks_predicted_predicted.RData")
    writeRData(networks, fileR)
    remove(fileR)
  }
  remove(netBc3, netC3, netAracne, netClr, netMrnet, netMrnetB)
  return(networks)
}