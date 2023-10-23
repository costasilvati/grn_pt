# Função para calcular métricas de comparação de matrizes
compareMatrices <- function(net, gold, toolName, threshold) {
  # Verifica se as matrizes têm as mesmas dimensões
  if (!identical(dim(net), dim(gold))) {
    stop("As matrizes net e gold devem ter as mesmas dimensões")
  }
  
  # Calcula os valores de FP, FN, TP e TN
  fp <- sum(net & !gold)
  fn <- sum(!net & gold)
  tp <- sum(net & gold)
  tn <- sum(!net & !gold)
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

