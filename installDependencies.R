installDependencies <- function(){
  install.packages("wrMisc")
  install.packages("Matrix")
  install.packages("c3net")
  install.packages("bc3net")
  install.packages("readr")
  install.packages("Matrix")
  install.packages("dplyr")
  install.packages("netdiffuseR")
  install.packages("tidyverse")
  install.packages("modelr")
  install.packages("tibble")
  install.packages("igraph")
  install.packages("stringr")
  install.packages("ennet")
  
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("minet")
  BiocManager::install("corto")
}