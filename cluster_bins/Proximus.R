##########
#Testing the Proximus algorithm for clustering discrete (binary) data
##########

library(cba)

#featureMatrix
featureMatrix <- read.table("/mnt/data0/noah/3C_methods/featureMatrix_idea031322/Granta519_0hr_SMC1-YY1-H3K27ac-H3K4me1-NIBPL_binned_res10000.matrix",
                            header=T,
                            sep="\t")
rownames(featureMatrix) <- featureMatrix$X
featureMatrix <- featureMatrix[,-c(1:4)]
featureMatrix <- as.matrix(mutate_all(featureMatrix, as.logical))

### cluster with proxiumus
# all default parameters
pr <- proximus(x=featureMatrix, max.radius = 1.75, debug=T)
approximations <- pr$a
