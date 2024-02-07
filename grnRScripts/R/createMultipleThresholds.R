# This fuction write RData to save memory
# Pattern name pathOut + _ + tool + threshold + filtred_geneGene.RData
createMultipleThresholds <- function(networks, 
                                      goldNet,
                                      tfList, 
                                      nSteps = 10,
                                      namesCol = c("fp","fn","tp","tn","accuracy","recall","precision","specificity","f_score","fdr", "AUPR", "AUC"),
                                     max_records = 10000, 
                                     writeData=FALSE, 
                                      pathOut = "."){
  names_net <- names(networks)
  results <- data.frame(matrix(0, 
                               nrow = (length(networks) * nSteps),
                               ncol = length(namesCol)
                               )
                        )
  i <- 1
  countRow <- 1
  #networksClean <- list()
  for (net in networks) {
    # -- remove connections < threshold
    minThreshold <- min(net)
    maxThreshold <- max(net)
    incrementValue = (maxThreshold- minThreshold)/nSteps
    nameList <- names_net[i]
    namesT <- NULL
    t <- 1
    cat(nameList, "- Threshold: ", minThreshold, " to ", maxThreshold, "increment: ", incrementValue)
    while (t <= nSteps) {
      cat("\n filter by threshold: ",minThreshold)
      netF <- edgeMatrixByTreshold(net,
                                   minThreshold,
                                   writeData = FALSE, 
                                   paste0(pathOut,
                                          "multipleThresholds/threshold/",
                                          nameList,"_", 
                                          minThreshold
                                          )
                                   )
      netClean <- removeGeneGeneConnection(tfList,
                                           netF,
                                           writeData, 
                                           paste0(pathOut,"multipleThresholds/geneGeneConnection/", 
                                                  nameList,
                                                  "_",
                                                  round(minThreshold, digits = 4)))
      #networksClean[[paste0(nameList,"_",t,"_",minThreshold)]] <- netClean
      results[countRow,] <- compareMatrices(netClean, 
                                            goldNet, 
                                            nameList, 
                                            minThreshold,
                                            maxRecords = max_records)
      minThreshold <- minThreshold + incrementValue
      countRow <- countRow + 1
      t <- t+1
    }
    i <- i+1
  }
  names(results) <- namesCol
  remove(names_net, minThreshold, nameList, namesT, i, t)
  #data <- list(results, networksClean)
  #return(data)
  return(results)
}