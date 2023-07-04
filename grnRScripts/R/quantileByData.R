quantileByData <- function(net, quantile = '100%'){
  netSorted <- sort(net)
  treshold <- quantile(netSorted)
  valueQuantile <- unname(treshold['100%']) # get value in named vector 
  return(valueQuantile)
}