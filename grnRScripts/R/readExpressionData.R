readExpressionData <- function(path){
  expData <- read.delim(path, sep = "\t", header = TRUE)
  return(expData)
}