# Função para calcular métricas de comparação de matrizes
compareMatrices <- function(net, gold, toolName, threshold,maxRecords) {
  # Verifica se as matrizes têm as mesmas dimensões
  if (!identical(dim(net), dim(gold))) {
    stop("As matrizes net e gold devem ter as mesmas dimensões")
  }
  
  # Limite o número de registros a serem comparados
  num_records <- min(maxRecords, nrow(net) * ncol(net))
  
  # Máscara para o triângulo superior
  upper_triangle <- upper.tri(net)
  
  # Aplicar a máscara e calcular TP, FP, TN e FN usando operações de soma bit a bit
  tp <- sum(net[upper_triangle] & gold[upper_triangle])
  fp <- sum(net[upper_triangle] & !gold[upper_triangle])
  tn <- sum(!net[upper_triangle] & !gold[upper_triangle])
  fn <- sum(!net[upper_triangle] & gold[upper_triangle])
  # Calcula as métricas de desempenho
  accuracy <- (tp + tn) / (tp + tn + fp + fn)
  recall <- tp / (tp + fn) # == TPR
  precision <- tp / (tp + fp)
  specificity <- tn / (tn + fn)
  f_score <- 2 * (precision * recall) / (precision + recall)
  fdr = fp/(fp + tn)
  # Cria um data.frame com os resultados
  results <- data.frame(
    tool = toolName,
    threshold = round(threshold, digits = 4),
    fp = fp,
    fn = fn,
    tp = tp,
    tn = tn,
    accuracy = accuracy,
    recall = recall,
    precision = precision,
    specificity = specificity,
    f_score = f_score,
    fdr = fdr
  )
  return(results)
}


