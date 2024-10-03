# Analyze the Three Toy Challenge Datasets
# Fitting DPMMs to learn the density
library(tidyverse)
library(grid)
library(gridExtra)
library(dirichletprocess)
library(salso)
library(mcclust.ext)
library(matrixStats)
library(kneedle)
library(ggplot2)
source("R/dirichletprocess_custom_init.R")
source("R/density_clusterer.R")
source("R/ne_parts_pair_counting.R")
source("R/salso_custom.R")
source("R/credible_bounds.R")
source("R/rearrange_labels.R")
source("R/choose_lambda.R")
Rcpp::sourceCpp("src/subpartiton_min_max.cpp")
dir.create('output/toy_challenge', recursive = TRUE, showWarnings = FALSE)


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

# This Loop Will Create data-frames for each data-set that contain a variety of 
# information we would like to plot.
nsims <- 1000 #number of MCMC runs.

for(d in 1:length(toy_datasets)){
  
  cat("Loading Dataset: ", names(toy_datasets)[d], "\n")
  dataset_name <- names(toy_datasets)[d]
  x <- scale(as.matrix(toy_datasets[[d]]))
  xx <- yy <- seq(-3,3, length.out = 30)
  xy_grid <- as.matrix(expand.grid(xx,yy))
  
  plot_obs <- toy_datasets[[d]]
  plot_grid <- data.frame(x = xy_grid[,1],
                          y = xy_grid[,2])
  
  g0Priors <- list(mu0 = rep(0,ncol(x)),
                   Lambda = diag(ncol(x)),
                   kappa0 = 0.1, #0.01,
                   nu = 10) # 200
  
  # Fit the DPMM model if not done so already.
  R.cache::evalWithMemoization({
    set.seed(12345)
    
    cat("Fitting DPMM Model\n")
    dp_mod <- DirichletProcessMvnormal(x, 
                                     numInitialClusters = 20,
                                     g0Priors = g0Priors,
                                     alphaPriors = c(1,0.01))
    dp_mod <- custom_init(dp_mod)
    dp_mod <- Fit(dp_mod, nsims)
  }, key=list(g0Priors=g0Priors, x=x, nsims=nsims))
  
  # Extract the Density Samples at the Observations and on a Grid
  fn_samps_obs <- matrix(nrow = nsims, ncol = nrow(x))
  for(s in 1:nsims)
    fn_samps_obs[s,] <- PosteriorFunction(dp_mod, s)(x)
  plot_obs$f_pe <- colMeans(fn_samps_obs)
  
  fn_samps_grid <- matrix(nrow = nsims, ncol = nrow(xy_grid))
  for(s in 1:nsims)
    fn_samps_grid[s,] <- PosteriorFunction(dp_mod, s)(xy_grid)
  plot_grid$f_pe <- colMeans(fn_samps_grid)
  rm(fn_samps_grid)
  
  cat("Fitting DPMM Model with SALSO\n")
  # Get the Model-Based Cluster Allocations and WG Vertical Credible Bounds
  mixture_clustering_samps <- matrix(nrow = nsims, ncol = nrow(x))
  for(s in 1:nsims)
    mixture_clustering_samps[s,] <- dp_mod$labelsChain[[s]]
  
  mixture_pe <- salso::salso(mixture_clustering_samps, 
                             maxZealousAttempts = 25)
  mixture_bounds <- model_based_credible_bounds(mixture_pe,
                                                mixture_clustering_samps)
  wg_mixture_bounds <- mcclust.ext::credibleball(mixture_pe, mixture_clustering_samps)
  plot_obs$mm_pe <- mixture_pe
  plot_obs$mm_vl <- mixture_bounds$lower
  plot_obs$mm_vu <- mixture_bounds$upper
  plot_obs$mm_vl_wg <- wg_mixture_bounds$c.lowervert[1,]
  plot_obs$mm_vu_wg <- wg_mixture_bounds$c.uppervert[1,]
  rm(dp_mod)
  
  # Determine the fraction of noise points using the elbow plot.
  # See the illustrations/gaussian_splitting_example.R for more details.
  Ef <- matrixStats::colMedians(fn_samps_obs)
  rank_and_den <- kneedle(rank(Ef), log(Ef), decreasing=FALSE, concave = TRUE, sensitivity = 2)
  #cut_quantile <- rank_and_den[1]/length(Ef)
  
  qplot(rank(Ef), log(Ef), xlab="ranks", ylab="sorted log(density)") + 
    geom_vline(xintercept=rank_and_den[1], color="red")
  ggsave(paste0("output/toy_challenge/level_selection_dpmm_", dataset_name, ".png"), 
         width = 10, height = 10)
  
  cat(sprintf("Calculating the cluster tree for %s\n", dataset_name))
  clusters <- level_set_clusters(x, Ef, cut_quantiles=seq(0,1,length.out=100))
  saveRDS(clusters, paste0("output/toy_challenge/clustree_dpmm_", dataset_name, ".rds"))
  saveRDS(Ef, paste0("output/toy_challenge/density_pe_dpmm_", dataset_name, ".rds"))
  
  #clustree(clusters, prefix="q")
  #ggsave(paste0("output/toy_challenge/clustree_dpmm_", dataset_name, ".png"), 
  #      width = 10, height = 10)
  
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
    plot_obs[[sprintf("db_pe_%.2f", q)]] <- res$pe
    plot_obs[[sprintf("db_vl_%.2f", q)]] <- res$lower_bound
    plot_obs[[sprintf("db_vu_%.2f", q)]] <- res$upper_bound
  }
  
  #plot_obs$db_pe <- res$pe
  #plot_obs$db_vl <- res$bounds$lower
  #plot_obs$db_vu <- res$bounds$upper
  #plot_obs$db_noise_frac <- quantiles[[dataset_name]]
  
  saveRDS(plot_obs, paste0("output/toy_challenge/plot_obs_dpmm_", dataset_name, ".rds"))
  saveRDS(plot_grid, paste0("output/toy_challenge/plot_grid_dpmm_", dataset_name, ".rds"))
  
  rm(fn_samps_obs)
}
