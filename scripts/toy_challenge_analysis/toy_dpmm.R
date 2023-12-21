# Analyze the Three Toy Challenge Datasets
# Fitting DPMMs to learn the density
library(tidyverse)
library(grid)
library(gridExtra)
library(dirichletprocess)
library(salso)
library(mcclust.ext)
library(matrixStats)
source("R/dirichletprocess_custom_init.R")
source("R/density_clusterer.R")
source("R/ne_parts_pair_counting.R")
source("R/salso_custom.R")
source("R/credible_bounds.R")
source("R/rearrange_labels.R")
source("R/choose_lambda.R")

Rcpp::sourceCpp("src/subpartiton_min_max.cpp")
dir.create('output/toy_challenge', recursive = TRUE, showWarnings = FALSE)

toy_datasets <- list(
    tsne = "tsne",
    two_moons = "two_moons",
    circles = "circles",
    tsne2 = "tsne"
)

cut_quantiles <- list(tsne =  0.2,
                     two_moons = 0.08,
                     circles = 0.0,
                     tsne2 = 0.1)

random_seed <- Sys.getenv("SLURM_ARRAY_TASK_ID")
selected_dataset <- 1 + random_seed %% 4


set.seed(random_seed)

two_moons <- readRDS("data/clean_data/two_moons.rds")
circles <- readRDS("data/clean_data/noisy_circles.rds")
tsne <- readRDS("data/clean_data/tsne.rds")

visualize_cluster_tree <- TRUE

# This Loop Will Create dataframes for each dataset that contain a variety of 
# information we would like to plot.
nsims <- 1000

cat("Loading Dataset: ", names(toy_datasets)[selected_dataset], "\n")
dataset_name <- names(toy_datasets)[selected_dataset]
data <- readRDS(sprintf("data/clean_data/%s", toy_datasets[[selected_dataset]]))

x <- scale(as.matrix(data))
xx <- yy <- seq(-3,3, length.out = 30)
xy_grid <- as.matrix(expand.grid(xx,yy))
  
plot_obs <- data
plot_grid <- data.frame(x = xy_grid[,1],
                        y = xy_grid[,2])
  
g0Priors <- list(mu0 = rep(0,ncol(x)),
                Lambda = diag(ncol(x)),
                kappa0 = 0.1, #0.01,
                nu = 10) # 200

# Fit the DPMM  Model
cat("Fitting DPMM Model\n")
dp_mod <- DirichletProcessMvnormal(x, 
                                  numInitialClusters = 20,
                                  g0Priors = g0Priors,
                                  alphaPriors = c(1,0.01))
dp_mod <- custom_init(dp_mod)
dp_mod <- Fit(dp_mod, nsims)
  
# Extract the Density Samples at the Observations and on a Grid
fn_samps_obs <- matrix(nrow = nsims, ncol = nrow(x))
for(s in 1:nsims)
  fn_samps_obs[s,] <- PosteriorFunction(dp_mod, s)(x)

plot_obs$f_pe <- colMedians(fn_samps_obs)
  
fn_samps_grid <- matrix(nrow = nsims, ncol = nrow(xy_grid))

for(s in 1:nsims)
  fn_samps_grid[s,] <- PosteriorFunction(dp_mod, s)(xy_grid)

plot_grid$f_pe <- colMedians(fn_samps_grid)
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
  
  
if (visualize_cluster_tree) {
  cat("Visualizing Cluster Tree for BALLET\n")
  Ef <- matrixStats::colMedians(fn_samps_obs)
  clusters <- level_set_clusters(x,Ef, cut_quantiles=seq(0,1,by=0.1))
  clustree(clusters, prefix="q")
  ggsave(paste0("output/toy_challenge/clustree_", dataset_name, ".png"), 
          width = 10, height = 10)
}
  
cat("Finding BALLET density based-clusters\n")
# Get the Density-Based Cluster Allocations and Our Credible Bounds
density_clustering_samps <- 
    density_based_clusterer(x, fn_samps_obs,
                            cut_quantile = cut_quantiles[[selected_dataset]])
pst <- compute_pst(density_clustering_samps)
pdt <- compute_pdt(density_clustering_samps)
density_pe <- salso_custom(density_clustering_samps, pst, pdt)
density_bounds <- credible_ball_bounds_active_inactive(x, density_pe, density_clustering_samps)
  
plot_obs$db_pe <- density_pe
plot_obs$db_vl <- density_bounds$lower
plot_obs$db_vu <- density_bounds$upper
  
  
saveRDS(plot_obs, paste0("output/toy_challenge/plot_obs_dpmm_", dataset_name, ".rds"))
saveRDS(plot_grid, paste0("output/toy_challenge/plot_grid_dpmm_", dataset_name, ".rds"))
