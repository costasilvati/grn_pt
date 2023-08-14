quantileByData <- function(net, quantil = '75%'){
  netSorted <- sort(net)
  treshold <- quantile(netSorted)
  valueQuantile <- unname(treshold[quantil]) # get value in named vector 
  return(valueQuantile)
}