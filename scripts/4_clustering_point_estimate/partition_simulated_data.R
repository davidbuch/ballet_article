library(tidyverse)
library(dirichletprocess)
source("R/density_clusterer.R")

simulated_datasets <- c("aniso", 
                        "blobs",
                        "two_moons")

for(ds_name in simulated_datasets){
  fsamps <- readRDS(sprintf("fitted_models/simulated_data/%s_dpmm_c1.rds", ds_name))
  Y <- dp_fit$data
  N <- nrow(Y)
  P <- ncol(Y)
  S <- length(dp_fit$alphaChain)
  
  burn <- floor(S/2)
  S <- S - burn
  
  density_clustering_samps <- matrix(nrow = S, ncol = N)
  mixture_clustering_samps <- matrix(nrow = S, ncol = N)
  for(s in 1:S){
    df <- PosteriorFunction(dp_fit, ind = burn + s)
    density_clustering_samps[s,] <- density_based_clusterer(Y, df)
    mixture_clustering_samps[s,] <- dp_fit$clusterLabels
  }
  
  clustering_samps %>%
    saveRDS(sprintf("fitted_models/simulated_data/%s_partitions_c1.rds", ds_name))
}


library(salso)
min_vi_part <- salso(clustering_samps)

ggplot() + 
  geom_point(aes(x = Y[,1], y = Y[,2], color = factor(min_vi_part)))

for(i in sample(1:S, 5)){
  p <- ggplot() + 
    geom_point(aes(x = Y[,1], y = Y[,2], color = factor(clustering_samps[i,]))) + 
    guides(color = "none")
  print(p)
}


