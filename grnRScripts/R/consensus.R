consensus <- function(networks, gold, namesCol, minConsenso = 3, writeData = TRUE, pathOut="."){
  dimensao_padrao <- dim(networks[[1]])
  for (i in 2:length(networks)) {
    if (!identical(dim(networks[[i]]), dimensao_padrao)) {
      stop("Todas as matrizes devem ter as mesmas dimensões.")
    }
  }
  toolNames<- names(networks)
  matriz_consenso <- mergeAdjacencyMatrices2(networks, minConsenso)
  resultsCons <- data.frame(matrix(0, 
                               nrow = (length(networks)),
                               ncol = length(namesCol)))
  cat("After consider ", minConsenso, "tools apointments. Number of Connections: ", sum(matriz_consenso == 1), "\n")
  resultsCons <- compareMatrices(matriz_consenso, gold, minConsenso, NaN)
  if(writeData){
    writeRData(matriz_consenso, paste0(pathOut, "netConsenso.RData"))
  }
  names(resultsCons) <- namesCol
  analisys <- list(net = matriz_consenso, result = resultsCons)
  remove(dimensao_padrao, i, matriz_consenso, resultsCons)
  return(analisys)
}