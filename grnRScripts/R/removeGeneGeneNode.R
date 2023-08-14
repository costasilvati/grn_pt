removeGeneGeneNode <- function(tfList, net, writeData = FALSE, pathOut = "."){
  message(paste((sum(net == 1)),"Node found without TF filter"))
  nomes_linhas <- row.names(net)
  nomes_colunas <- colnames(net)
  
  for (linha in nomes_linhas) {
    for (coluna in nomes_colunas) {
      if (net[linha, coluna] == 1 && (linha %in% tfList || coluna %in% tfList)) {
        paste("A string", string_alvo, "ocorre na posição", linha, "-", coluna)
      }
    }
  }
  
  if(writeData){
    message("Writting data...")
    fileName <- paste0(pathOut,"filtred_geneGene.csv")
    writeNetworkCsv(net, fileName)
    fileR = paste0(pathOut,"filtred_geneGene.RData")
    writeRData(net, fileR)
    remove(fileR)
  }
  remove(nomes_colunas)
  remove(nomes_linhas)
  remove(fileName)
  return(net)
}