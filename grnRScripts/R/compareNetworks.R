compareNetworks <- function(gold, net){
  gold_matrix <- as.matrix(gold)
  net_matrix <- as.matrix(net)
  TP <- sum(gold_matrix & net_matrix)
  TN <- sum(!gold_matrix & !net_matrix)
  FP <- sum(!gold_matrix & net_matrix)
  FN <- sum(gold_matrix & !net_matrix)
    
  confusion_matrix <- matrix(c(TP, FN, FP, TN), nrow = 2, byrow = TRUE,
                             dimnames = list(c("Positive", "Negative"), c("TRUE", "FALSE")))
  remove(gold_matrix)
  remove(net_matrix)
  return(confusion_matrix)
}

# Função para calcular matriz de confusão
conf_mat <- function(gold, net) {
  tp <- sum(gold * net)
  fp <- sum(net) - tp
  fn <- sum(gold) - tp
  tn <- sum((gold == 0) & (net == 0))
  return(matrix(c(tp, fp, fn, tn), nrow = 2, byrow = TRUE))
}

# Função para calcular taxas de falso positivo e true positivo
fp_rate <- function(conf_mat) {
  return(conf_mat[2, 1] / sum(conf_mat[2, ]))
}
tp_rate <- function(conf_mat) {
  return(conf_mat[1, 1] / sum(conf_mat[1, ]))
}

# Função para calcular acurácia
acc <- function(conf_mat) {
  return((conf_mat[1, 1] + conf_mat[2, 2]) / sum(conf_mat))
}

# Função para gerar curva ROC
roc_curve <- function(gold, net) {
  conf_mat <- conf_mat(gold, net)
  fp_rate <- fp_rate(conf_mat)
  tp_rate <- tp_rate(conf_mat)
  return(list(fp_rate = fp_rate, tp_rate = tp_rate))
}

