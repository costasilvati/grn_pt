# This fuction write RData to save memory
# Pattern name pathOut + _ + tool + threshold + filtred_geneGene.RData
createMultipleThresholds <- function(networks, tfList, goldStandard, writeData=FALSE, pathOut = "."){
  names_net <- names(networks)
  netThresolds <- list()
  netThresolds <- names(names_net)
  i <- 1
  for (net in networks) {
    minThreshold <- min(net)
    nameList <- names_net[i]
    namesT <- NULL
    t <- 1
    while (minThreshold <= 1.0) {
      message(paste(nameList," filter by threshold: ",minThreshold))
      netF <- edgeMatrixByTreshold(net,minThreshold,writeData = TRUE, paste0(pathOut,nameList,"_"))
      message(paste(nameList," removing gene-gene connections"))
      removeGeneGeneConnection(tfList,netF,TRUE, paste0(pathOut,nameList,"_",minThreshold))
      minThreshold <- minThreshold + 0.1
      t <- t+1
    }
    i <- i+1
  }
  remove(names_net, netThresolds, minThreshold, nameList, namesT, i, t)
}