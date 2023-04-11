if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("minet")
library(minet)
# ARACNE
data(syn.data)
mim <- build.mim(syn.data,estimator="spearman")
net <- aracne(mim)

#MRNetB
mr_net <- mrnetb(mim)
