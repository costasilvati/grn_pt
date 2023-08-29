edgeMatrixByTreshold <- function(network, thresholdNet = 0.0, writeData = FALSE, pathOut = "."){
  net_df <- as.data.frame(network)
  colunas <- colnames(net_df)
  suppressWarnings(
    net2 <- net_df %>% 
    mutate_at(
      vars(one_of(colunas)),
      funs(case_when(
        . > thresholdNet ~ 1,
        . <= thresholdNet ~ 0)))
  )
  if(writeData){
    #fileName <- paste0(pathOut,"aracne_predicted_original.csv")
    #writeNetworkCsv(net, fileName)
    fileR = paste0(pathOut,"threshold.RData")
    writeRData(net2, fileR)
  }
  message(paste((sum(net2 == 1)),"Node found."))
  remove(net_df)
  remove(colunas)
  return(net2)
}
