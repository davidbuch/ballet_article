# Analyze the Three Toy Challenge Datasets
# Fitting Nearest Neighbor Dirichlet Mixtures to learn the density
library(plyr)  # preemptively load plyr before NNDP loads dplyr/plyr out of order
library(tidyverse)
library(NNDM) # devtools::install_github("shounakchattopadhyay/NN-DM")
source("R/density_clusterer.R")
source("R/ne_parts_pair_counting.R")
source("R/salso_custom.R")
source("R/credible_bounds.R")
Rcpp::sourceCpp("src/subpartiton_min_max.cpp")
dir.create('output/toy_challenge', recursive = TRUE, showWarnings = FALSE)

set.seed(1234)

two_moons <- readRDS("data/clean_data/two_moons.rds")
circles <- readRDS("data/clean_data/noisy_circles.rds")
tsne <- readRDS("data/clean_data/tsne.rds")

toy_datasets <- list(
  two_moons = two_moons,
  circles = circles,
  tsne = tsne
)

# This Loop Will Create dataframes for each dataset that contain a variety of 
# information we would like to plot.
nsims <- 1000
for(d in 1:length(toy_datasets)){
  dataset_name <- names(toy_datasets)[d]
  if(!file.exists(paste0("output/toy_challenge/plot_obs_nndm_", dataset_name, ".rds")) ||
     !file.exists(paste0("output/toy_challenge/plot_grid_nndm_", dataset_name, ".rds"))){
    x <- scale(as.matrix(toy_datasets[[d]]))
    xx <- yy <- seq(-3,3, length.out = 30)
    xy_grid <- as.matrix(expand.grid(xx,yy))
    
    plot_obs <- toy_datasets[[d]]
    plot_grid <- data.frame(x = xy_grid[,1],
                            y = xy_grid[,2])
    
    # Fit NNDM model and Extract Density Samples at the Observations and on a Grid
    p <- ncol(x)
    mu0 <- rep(0,p)
    nu0 <- 0.001
    gamma0 <- p
    
    res <- mNNDM(x = x, MC = nsims, k = 10, inputpt = x, 
                 mu0 = mu0, nu0 = nu0, gamma0 = p)
    fn_samps_obs <- t(res$f_stor)
    plot_obs$f_pe <- colMeans(fn_samps_obs)
    
    res <- mNNDM(x = x, MC = nsims, k = 10, inputpt = xy_grid, 
                 mu0 = mu0, nu0 = nu0, gamma0 = p)
    fn_samps_grid <- t(res$f_stor)
    plot_grid$f_pe <- colMeans(fn_samps_grid)
    plot_grid$f_s1 <- fn_samps_grid[round(nsims/2),]
    plot_grid$f_s2 <- fn_samps_grid[nsims,]
    rm(fn_samps_grid)
    
    # Get the Density-Based Cluster Allocations and Our Credible Bounds
    density_clustering_samps <- 
      density_based_clusterer(x, fn_samps_obs, 
                              cut_quantile = 0.125, 
                              split_err_prob = 0.01)
    rm(fn_samps_obs)
    pst <- compute_pst(density_clustering_samps)
    pdt <- compute_pdt(density_clustering_samps)
    density_pe <- salso_custom(density_clustering_samps, pst, pdt)
    density_bounds <- credible_ball_bounds_active_inactive(x, density_pe, density_clustering_samps)
    
    plot_obs$db_pe <- density_pe
    plot_obs$db_vl <- density_bounds$lower
    plot_obs$db_vu <- density_bounds$upper
    
    saveRDS(plot_obs, paste0("output/toy_challenge/plot_obs_nndm_", dataset_name, ".rds"))
    saveRDS(plot_grid, paste0("output/toy_challenge/plot_grid_nndm_", dataset_name, ".rds"))
  }
}
