# Analyze the Three Toy Challenge Datasets
# Fitting Adaptive Polya Trees to learn the density
library(tidyverse)
library(PTT) # devtools::install_github("MaStatLab/PTT")
library(matrixStats)
library(kneedle)
library(ggplot2)
source("R/scale_box.R")
source("R/fn_from_apt.R")
source("R/density_clusterer.R")
source("R/ne_parts_pair_counting.R")
source("R/salso_custom.R")
source("R/credible_bounds.R")
source("R/choose_lambda.R")
Rcpp::sourceCpp("src/subpartiton_min_max.cpp")
dir.create('output/toy_challenge', recursive = TRUE, showWarnings = FALSE)

set.seed(1234)

two_moons <- readRDS("data/clean_data/two_moons.rds")
circles <- readRDS("data/clean_data/noisy_circles.rds")
tsne <- readRDS("data/clean_data/tsne.rds")

toy_datasets <- list(
  tsne = tsne,
  two_moons = two_moons,
  circles = circles
)

# An NA means: determine the level using the elbow heuristic.
quantiles <- list(tsne =  c(0.15, 0.10, 0.05, NA),
                  two_moons = c(0.10, 0.05, NA),
                  circles = c(0.10, 0.05, NA))

# This Loop Will Create dataframes for each dataset that contain a variety of 
# information we would like to plot.
nsims <- 1000
for(d in 1:length(toy_datasets)){
  cat("Loading Dataset: ", names(toy_datasets)[d], "\n")
  dataset_name <- names(toy_datasets)[d]
  x <- scale_box(toy_datasets[[d]])
  xx <- yy <- seq(0.00001,0.99999, length.out = 30)
  xy_grid <- as.matrix(expand.grid(xx,yy))
  
  plot_obs <- toy_datasets[[d]]
  
  plot_grid <- data.frame(x = xy_grid[,1], y = xy_grid[,2])
  
  # Fit the APT Model and Extract Density samples for both the observed x and the grid
  set.seed(1234)
    
  cat("Fitting APT Model\n")
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
  
  # Determine the fraction of noise points using the elbow heuristic.
  # See the illustrations/gaussian_splitting_example.R for more details.
  # Finding the elbow here needed a higher value of "sensitivity".
  
  Ef <- matrixStats::colMedians(fn_samps_obs)
  rank_and_den <- kneedle(rank(Ef), log(Ef), decreasing=FALSE, concave = TRUE)
  cut_quantile <- rank_and_den[1]/length(Ef)
  
  qplot(rank(Ef), log(Ef), xlab="ranks", ylab="sorted log(density)") + 
    geom_vline(xintercept=rank_and_den[1], color="red")
  ggsave(paste0("output/toy_challenge/level_selection_apt_", dataset_name, ".png"), 
         width = 10, height = 10)
  
  cat(sprintf("Calculating the cluster tree for %s\n", dataset_name))
  clusters <- level_set_clusters(x, Ef, cut_quantiles=seq(0,1,length.out=100))
  saveRDS(clusters, paste0("output/toy_challenge/clustree_apt_", dataset_name, ".rds"))
  saveRDS(Ef, paste0("output/toy_challenge/density_pe_apt_", dataset_name, ".rds"))
  
  run_ballet <- function(cut_quantile) {
    
    cat(sprintf("Finding BALLET density-based clusters 
                with %.2f%% noise.\n", 100*cut_quantile))
    # Get the Density-Based Cluster Allocations and Our Credible Bounds
    density_clustering_samps <- 
      density_based_clusterer(x, fn_samps_obs,
                              cut_quantile = cut_quantile)
    pst <- compute_pst(density_clustering_samps)
    pdt <- compute_pdt(density_clustering_samps)
    density_pe <- salso_custom(density_clustering_samps, pst, pdt)
    density_bounds <- credible_ball_bounds_active_inactive(x, density_pe, 
                                                           density_clustering_samps)
    
    lower_bound=density_bounds$lower
    upper_bound=density_bounds$upper
    
    attr(density_pe, "noise_frac") <- cut_quantile
    attr(lower_bound, "noise_frac") <- cut_quantile
    attr(upper_bound, "noise_frac") <- cut_quantile
    
    list(pe=density_pe,
         lower_bound=lower_bound,
         upper_bound=upper_bound)
    
  }
  
  # Run BALLET with different choices of level
  for(i in 1:length(quantiles[[dataset_name]])) {
    q <- quantiles[[dataset_name]][i]
    if(is.na(q)) {
      # An NA value means use our heuristic threshold.
      q <- rank_and_den[1]/length(Ef)
    }
    res <- run_ballet(q) 
    
    # Store the quantile value in the column name..  
    plot_obs[[sprintf("db_pe_%.4f", q)]] <- res$pe
    plot_obs[[sprintf("db_vl_%.4f", q)]] <- res$lower_bound
    plot_obs[[sprintf("db_vu_%.4f", q)]] <- res$upper_bound
  }
  
  saveRDS(plot_obs, paste0("output/toy_challenge/plot_obs_apt_", dataset_name, ".rds"))
  saveRDS(plot_grid, paste0("output/toy_challenge/plot_grid_apt_", dataset_name, ".rds"))
  
  rm(fn_samps_obs)
}
