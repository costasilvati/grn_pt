readExpressionData <- function(path, demlimiterChar = "\t"){
  expData <- read.delim(path, sep =  demlimiterChar, header = TRUE)
  return(expData)
}