# Read a list nodes and return a adjacenc matrix
listToMatrix <- function(pathListFile, 
                         colRowNames,
                         headerExsit=FALSE, 
                         separator=",", 
                         dimension=4511, 
                         writeRData= FALSE, 
                         pathOut = ".") {
  listData <- read.csv(pathListFile, header = headerExsit, sep = separator) # edge list

  matrixData <- matrix(data = 0, nrow = dimension, ncol = dimension) # matrix
  row.names(matrixData) <- colRowNames
  colnames(matrixData) <- colRowNames
  node1 <- listData[1]
  node2 <- listData[2]
  edge <- listData[3]
  line <- 1
  for (n1 in node1[,1]) {
    matrixData[n1,node2[line, 1]] <- edge[line, 1]
    line <- line +1
  }
  if(writeRData){
    writeRData(listData, paste0(pathOut,"_List.RData"))
    writeRData(matrixData, paste0(pathOut, "_Matrix.RData"))
  }
  sum(matrixData == 1)
  sum(listData$V3 == 1)
  remove(line, node1, node2, edge, n1)
  return(matrixData)
}