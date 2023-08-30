inferenceNetwork <- function(expressionData, writeData = FALSE, pathOut = "."){
  # -- gera matriz transposta
  expressionDataTransposto <- t(expressionData)
  colnames(expressionDataTransposto) <- rep(1:length(row.names(expressionData)))
  #--- discretizar para métodos MI  # default do infotheo função discretize()
  message("Executing unsupervised data discretization...")
  expressionTranspDiscret <- infotheo::discretize(expressionDataTransposto)
  message("Executing NetBc3...")
  netBc3 <- runBc3Net(expressionDataTransposto,writeData, pathOut)
  message("Executing NetC3...")
  netC3 <- runC3net(expressionTranspDiscret, writeData, pathOut)
  message("Executing ARACNE...")
  netAracne <- runAracne(expressionTranspDiscret, writeData, pathOut)
  message("Executing CLR...")
  netClr <- runClr(expressionTranspDiscret, writeData, pathOut)
  message("Executing MRNet...")
  netMrnet <- runMRnet(expressionTranspDiscret, writeData, pathOut)
  message("Executing MRNetB...")
  netMrnetB <- runMrnetB(expressionTranspDiscret, writeData, pathOut)
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