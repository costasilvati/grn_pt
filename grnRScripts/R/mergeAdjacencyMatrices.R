mergeAdjacencyMatrices <- function(networks, minConsenso = 3){
  n <- length(networks)
  if (n == 0) {
    stop("A lista de matrizes de adjacência está vazia.")
  }
  matrix_size <- dim(networks[[1]])
  result_matrix <- matrix(0, nrow = matrix_size[1], ncol = matrix_size[2])
  row.names(result_matrix) <- row.names(networks[[1]])
  colnames(result_matrix) <- colnames(networks[[1]])
  for (i in 1:matrix_size[1]) {
    for (j in 1:matrix_size[2]) {
      count <- sum(sapply(networks, function(x) x[i, j]))
      if (count >= minConsenso) {
        result_matrix[i, j] <- 1
      }
    }
  }
  return(result_matrix)
}
#-- use a 2
mergeAdjacencyMatrices2 <- function(networks, minConsenso = 3) {
  n <- length(networks)
  if (n == 0) {
    stop("A lista de matrizes de adjacência está vazia.")
  }
  
  matrix_size <- dim(networks[[1]])
  result_matrix <- matrix(0, nrow = matrix_size[1], ncol = matrix_size[2])
  row.names(result_matrix) <- row.names(networks[[1]])
  colnames(result_matrix) <- colnames(networks[[1]])
  
  # Soma de todas as matrizes na lista
  total_sum <- Reduce(`+`, networks)
  
  # Aplicar condição de minConsenso
  result_matrix[total_sum >= minConsenso] <- 1
  
  return(result_matrix)
}
