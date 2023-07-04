library(tidyverse)
library(tibble)
library(dplyr)
#' Recebe uma matriz de adjacencia 
#'
edgeListByTreshold <- function(net, treshold = 0.0){
  net_tb <- as_tibble(net)
  nf <- 0
  l<- 1
  genes <- colnames(net_tb)
  edgeList <- tibble(node1=rep("node",1), node2=rep("node",1), value=rep(0.0,1))
  for (c in 1:length(net_tb[1,])) {
    while (l <= length(genes)) {
      if(net_tb[c,l] >= treshold){
        edgeList <- edgeList %>% add_row(tibble::tibble_row(node1=genes[c], node2=genes[l], value= as.double(net_tb[c,l])))
        nf <- nf + 1
      }
      l <- l + 1
    }
  }
  edgeList <- edgeList %>% filter(!dplyr::row_number() %in% c(1))
  message(paste(nf," Nodes found"))
  return(edgeList)
}
