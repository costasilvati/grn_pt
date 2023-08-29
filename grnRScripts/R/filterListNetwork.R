filterListNetwork <- function(networks, writeData = FALSE, pathOut = "."){
  names_net <- names(networks)
  networks_filtred <- list(names_net)
  cont <- 1
  for (net in networks) {
    message(names_net[cont])
    thresholdNetwork <- quantileByData(net)
    networks_filtred[[cont]] <- edgeMatrixByTreshold(network=net, 
                                                     threshold= thresholdNetwork, 
                                                     writeData = writeData, 
                                                     pathOut = paste(pathOut, names_net[cont])
    )
    cont <- cont + 1
  }
  names(networks_filtred) <- names_net 
  return(networks_filtred)
}