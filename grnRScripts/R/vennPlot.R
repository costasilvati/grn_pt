
encontra_uns <- function(matriz, gold) {
  posicoes <- which(matriz == 1, arr.ind = TRUE)
  posicoes <- posicoes[posicoes[,1] <= posicoes[,2], ]
  conName <- character(dim(posicoes)[1])
  rows <- row.names(matriz)
  cols <- colnames(matriz)
  tp <- 1
  for(i in 1:dim(posicoes)[1]){
    if(gold[posicoes[i,1], posicoes[i,2]] == 1){
      conName[tp] <- paste0(rows[posicoes[i,1]], "-",cols[posicoes[i,2]])
      tp <- tp + 1
    }else{
      # cat(rows[posicoes[i,1]], "-",cols[posicoes[i,2]], "FP \n")
    }
  }
  row.names(posicoes) <- conName
  return(conName)
}

makeUpsetPlot <- function(netoworksSet, 
                          condition = "condition_default",
                          goldStndardNetwork, 
                          pathOut = ".",
                          writeData = TRUE){
  resultados <- list()
  for(i in 1:length(networksSet)) {
    posicoes_uns <- encontra_uns(networksSet[[i]],goldStndardNetwork)
    resultados[[i]] <- posicoes_uns
  }
  names(resultados) <- names(networksSet)
  todas_posicoes <- unique(unlist(resultados))
  df <- data.frame(matrix(ncol = length(resultados), 
                          nrow = length(todas_posicoes)))
  rownames(df) <- todas_posicoes
  colnames(df) <- names(resultados)
  for(i in 1:length(resultados)) {
    df[resultados[[i]], i] <- 1
  }
  df[is.na(df)] <- 0
  linhas_G <- grep("^G", rownames(df))
  df2 <- subset(df, rownames(df) %in% rownames(df)[linhas_G])
  if(writeData){
    meu_grafico <- UpSetR::upset(df2, 
                                 sets = colnames(df2), 
                                 sets.bar.color = "#56B4E9",
                                 order.by = "freq", 
                                 empty.intersections = "on")
    pdf(paste0(pathOut,
               "/intersection_",
               length(names(resultados)),
               "_tools",condition ,
               ".pdf"), 
        width = 10, 
        height = 7.5, 
        onefile = FALSE)
    print(meu_grafico)
    dev.off()
    df2$`number of methods` <- rowSums(df2) 
    writeNetworkCsv(df2, 
                    paste0(pathOut,
                           "/intersection_",
                           length(names(resultados)),
                           "_tools",
                           condition ,
                           ".csv"))
  }else{
    df2$`number of methods` <- rowSums(df2) 
    UpSetR::upset(df2, 
                  sets = colnames(df2), 
                  sets.bar.color = "#56B4E9",
                  order.by = "freq", 
                  empty.intersections = "on", 
                  main = paste("Intersection",
                                length(names(resultados)),
                                " methods",condition))
  }
  return(df2)
}

