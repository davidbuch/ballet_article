# Analyze the Three Toy Challenge Datasets
# Fitting Adaptive Polya Trees to learn the density
library(tidyverse)
library(PTT) # devtools::install_github("MaStatLab/PTT")
source("R/scale_box.R")
source("R/fn_from_apt.R")
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

ballet_params <- list(
  two_moons = list(minPts = 5, cut_quantile = 0.08, split_err_prob = 0.01),
  circles = list(minPts = 5, cut_quantile = 0.025, split_err_prob = 0.01),
  tsne = list(minPts = floor((nrow(tsne)/nrow(circles))*5), 
              cut_quantile = 0.08, split_err_prob = 0.01)
)

# This Loop Will Create dataframes for each dataset that contain a variety of 
# information we would like to plot.
nsims <- 1000
for(d in 1:length(toy_datasets)){
  dataset_name <- names(toy_datasets)[d]
  x <- scale_box(toy_datasets[[d]])
  xx <- yy <- seq(0.00001,0.99999, length.out = 30)
  xy_grid <- as.matrix(expand.grid(xx,yy))
  
  plot_obs <- toy_datasets[[d]]
  
  plot_grid <- data.frame(x = xy_grid[,1], y = xy_grid[,2])
  
  # Fit the APT Model and Extract Density samples for both the observed x and the grid
  max.resol = 10
  res_obs <- apt(X = x, Xpred = x, n.post.samples = nsims, max.resol = max.resol)
  fn_samps_obs <- fn_from_apt(res_obs, x, max.resol) 
  fn_samps_obs <- fn_samps_obs / prod(attr(x, "scaled:range"))
  plot_obs$f_pe <- colMeans(fn_samps_obs)
  
  res_grid <- apt(X = x, Xpred = xy_grid, n.post.samples = nsims, max.resol = max.resol)
  fn_samps_grid <- fn_from_apt(res_grid, xy_grid, max.resol)
  
  # Transform the plot grid and density back to the original coordinates
  fn_samps_grid <- fn_samps_grid / prod(attr(x, "scaled:range"))
  plot_grid <- sweep(plot_grid, 2L, attr(x, "scaled:range"), '*', check.margin = FALSE)
  plot_grid <- sweep(plot_grid, 2L, attr(x, "scaled:lower"), '+', check.margin = FALSE)
  plot_grid$f_pe <- colMeans(fn_samps_grid)
  plot_grid$f_s1 <- fn_samps_grid[round(nsims/2),]
  plot_grid$f_s2 <- fn_samps_grid[nsims,]
  rm(fn_samps_grid)
  
  # Get the Density-Based Cluster Allocations and Our Credible Bounds
  density_clustering_samps <- 
    density_based_clusterer(x, fn_samps_obs,
                            minPts = ballet_params[[d]]['minPts'],
                            cut_quantile = ballet_params[[d]]['cut_quantile'], 
                            split_err_prob = ballet_params[[d]]['split_err_prob'])
  rm(fn_samps_obs)
  pst <- compute_pst(density_clustering_samps)
  pdt <- compute_pdt(density_clustering_samps)
  density_pe <- salso_custom(density_clustering_samps, pst, pdt)
  density_bounds <- credible_ball_bounds_active_inactive(x, density_pe, density_clustering_samps)
  
  plot_obs$db_pe <- density_pe
  plot_obs$db_vl <- density_bounds$lower
  plot_obs$db_vu <- density_bounds$upper
  
  saveRDS(plot_obs, paste0("output/toy_challenge/plot_obs_apt_", dataset_name, ".rds"))
  saveRDS(plot_grid, paste0("output/toy_challenge/plot_grid_apt_", dataset_name, ".rds"))
  
}
