compareNetworks <- function(gEdges, cEdges ){
  if(length(gEdges) == 3 | length(cEdges) == 3){
    for (l in 1:nrow(cEdges)) {
      for (gl in 1:nrow(gEdges)) {
        if(gEdges[l,3] == 1){
          node12 = cEdges[l,1] == gEdges[gl,1] & cEdges[l,2] == gEdges[gl,2] 
          node21 = cEdges[l,1] == gEdges[gl,2] & cEdges[l,2] == gEdges[gl,1] 
          if( node12 | node21 ){
            #node <- paste0(gEdges[gl,1], paste0(" - ", gEdges[gl,2]))
            paste0(gEdges[gl,1], paste0(" - ", gEdges[gl,2]))
          } 
        }
      }
    }
  }else{
    message("Gold Standard and other Network data, need 3 columns")
    return(NULL)
  }
}