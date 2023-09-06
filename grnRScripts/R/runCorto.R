runCorto <- function(dataExpression, tfList, writeCsv=FALSE, pathOut="."){
  
}
# obrigatório ter tfList
load(system.file("extdata","inmat.rda",package="corto"))
inmat[1:5,1:5]
dim(inmat)
load(system.file("extdata","centroids.rda",package="corto"))
centroids[15]
regulon<-corto(expressionDataTransposto,centroids=c(tfList$X1),nbootstraps=10,p=1e-30,nthreads=2)

remove(inmat, centroids, regulon)