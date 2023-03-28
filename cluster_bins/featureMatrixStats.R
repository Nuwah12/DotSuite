#1: Independence of features
# Test covariance between all features:
covs <- c()
for(i in seq(1:5)){
  for(j in seq(1:5)){
    covs <- c(covs, cov(featureMatrix[,i], featureMatrix[,j]))
  }
}
