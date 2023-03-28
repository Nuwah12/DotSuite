##########  
# Hierarchical clustering methods
# 3/27/23 - having trouble ingesting entiure matrix into agnes() clustering function
##########

library(tidyverse)
library(cluster)
library(factoextra)
library(dendextend)

featureMatrix <- read.table("/mnt/data0/noah/3C_methods/featureMatrix_idea031322/Granta519_0hr_SMC1-YY1-H3K27ac-H3K4me1-NIBPL_binned_res10000.matrix",
                            header=T,
                            sep="\t")
rownames(featureMatrix) <- featureMatrix$X
featureMatrix <- featureMatrix[,-c(1:4)]


##### agglomerative
# Jaccard distance
agg.complete <- dist(featureMatrix, method = "binary")
