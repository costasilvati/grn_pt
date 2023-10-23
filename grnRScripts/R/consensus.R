consensus <- function(networks, gold, namesCol, minConsenso = 3, writeData = TRUE, pathOut="."){
  dimensao_padrao <- dim(networks[[1]])
  for (i in 2:num_matrizes) {
    if (!identical(dim(networks[[i]]), dimensao_padrao)) {
      stop("Todas as matrizes devem ter as mesmas dimensões.")
    }
  }
  toolNames<- names(networks)
  matriz_consenso <- mergeAdjacencyMatrices(networks, minConcenso)
  cat("After consensus: Number of Connections: ", sum(matriz_consenso == 1), "\n")
  resultsCons <- data.frame(matrix(0, 
                               nrow = (length(networks)),
                               ncol = length(namesCol)))
  cat("After ", toolNames[i], "Number of Connections: ", sum(matriz_consenso == 1), "\n")
  resultsCons <- compareMatrices(matriz_consenso, gold, toolNames[i], NaN)
  
  if(writeData){
    writeRData(matriz_consenso, paste0(pathOut, "netConsenso.RData"))
  }
  names(resultsCons) <- namesCol
  analisys <- list(net = matriz_consenso, result = resultsCons)
  remove(dimensao_padrao, i, num_matrizes, matriz_consenso, resultsCons)
  return(analisys)
}