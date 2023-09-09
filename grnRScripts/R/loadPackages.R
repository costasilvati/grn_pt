loadPackages <- function(){
  dependenciesPakcage = c("corto", 
                           "stringr",
                           "bc3net",
                           "c3net",
                           "wrMisc",
                           "Matrix",
                           "readr",
                           "dplyr",
                           "netdiffuseR",
                           "tidyverse",
                           "minet",
                           "tibble",
                           "igraph",
                          "progress"
                          #,"ennet"
                          )
  invisible(lapply(dependenciesPakcage , library, character.only = TRUE))
  remove(dependenciesPakcage)
}