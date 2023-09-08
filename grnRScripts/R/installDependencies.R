installDependencies <- function(){
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BioconductorPackages = c("minet", "corto")
  
  for (biocPck in BioconductorPackages) {
    if (!requireNamespace(biocPck, quietly = TRUE)) {
      BiocManager::install(biocPck)
    }
  }
  
  CranPackages <- c("bc3net", 
               "c3net", 
               "wrMisc", 
               "ennet", 
               "netdiffuseR", 
               "modelr", 
               "tidyverse", 
               "tibble", 
               "igraph", 
               "stringr")
  
  # Verifique se cada pacote está instalado e instale-o se não estiver
  for (cranPck in CranPackages) {
    if (!requireNamespace(cranPck, quietly = TRUE)) {
      install.packages(cranPck)
    }
  }
}