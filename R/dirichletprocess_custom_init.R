custom_init <- function(dpObj){
  # Need to update
  # - clusterLabels
  # - pointsPerCluster
  # - clusterParameters
  
  # Fit kmeans
  K <- dpObj$numberClusters
  Y <- dpObj$data
  km_fit <- kmeans(x = Y, centers = K, nstart = 5)
  
  # Update clusterLabels
  clus_labs <- km_fit$cluster
  dpObj$clusterLabels <- clus_labs
  
  # Update pointsPerCluster
  dpObj$pointsPerCluster <- km_fit$size
  
  # Update clusterParameters (mu and sig)
  if(ncol(Y) > 1){
    for(k in 1:K){
      Y_k <- Y[clus_labs == k, , drop = FALSE]
      
      dpObj$clusterParameters$mu[,,k] <- colMeans(Y_k)
      
      # Only update the covariance if it is nonsingular
      sig <- cov(Y_k)
      if( !any(is.na(sig)) && !(det(sig) == 0)){
        dpObj$clusterParameters$sig[,,k] <- sig
      }
    }
  }else{
    for(k in 1:K){
      Y_k <- Y[clus_labs == k, , drop = FALSE]
      
      dpObj$clusterParameters$mu[k] <- colMeans(Y_k)
      
      # Only update the covariance if it is nonsingular
      sig <- cov(Y_k)
      if( !any(is.na(sig)) && !(det(sig) == 0)){
        dpObj$clusterParameters$sig[,,k] <- sig
      }
    }    
  }
  
  return(dpObj)
}