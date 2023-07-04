source("writeRData.R")
filterNetworkByTreshold <- function(net,
                                    writeRData = FALSE, 
                                    writeCsv = FALSE,
                                    pathOut = "../data/",
                                    tool = "na_"){
  treshold <- quantileByData(net)
  print(treshold)
  for (coluna in 1:length(net[1,])) {
    countRemoved <- net[coluna,net[coluna,] < treshold] <- 0
    countNode <- net[coluna,net[coluna,] >= treshold] <- 1
  }
  if(writeCsv){
    writeNetworkCsv(net, paste0(pathOut,paste0(tool,"_filtred.csv")))
  }
  if(writeRData){
    fileR = paste0(pathOut)
    writeRData(net, paste0(pathOut,paste0(tool,"_filtred.RData")))
  }
  return(net)
}