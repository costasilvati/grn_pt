# Os dados devem vir transpostos e com os nomes das colunas como nomes das linhas, por ex:
# >data <-read.csv('/mnt/dados/murilo/Network3_expression_data_T.csv')
# >rownames(data) <- data[, 1]
# >data <- data[, -1]

# Os tfs sao obrigatorios e devem estar como characters, por ex:
# >tfs <- read.delim('/mnt/dados/murilo/Network3_transcription_factors.tsv', stringsAsFactors = FALSE)
# >tfs <- as.character(unlist(tfs))

runCorto <- function(dataExpression, tfList, writeCsv=FALSE, pathOut="."){
  regulon<-corto(dataExpression,centroids=c(tfList$X1),nbootstraps=10,p=1e-30,nthreads=2)
  net <- create_and_populate_matrix(dataExpression,regulon)
  if(writeCsv){
    fileR = paste0(pathOut,"CORTO_predicted_original.RData")
    #save(net, file = fileR)
    writeRData(net, fileR)
  }
  return(net)
}

create_and_populate_matrix <- function(dataExpression,regulon){
  nomes <- rownames(dataExpression)
  num_nomes <- length(nomes)
  matriz_adjacencia <- matrix(0, nrow = num_nomes, ncol =num_nomes)
  rownames(matriz_adjacencia) <- nomes
  colnames(matriz_adjacencia) <- nomes
 
  for( tf in names(regulon)){

    for( gene in names(regulon[tf][[1]]$tfmode)){
      matriz_adjacencia[tf,gene] <- regulon[tf][[1]]$tfmode[[gene]] 
    }

  }
  return(matriz_adjacencia)
}
