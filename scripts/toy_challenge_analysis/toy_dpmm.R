# Analyze the Three Toy Challenge Datasets
# Fitting DPMMs to learn the density
library(tidyverse)
library(grid)
library(gridExtra)
library(concaveman)
library(dirichletprocess)
library(salso)
library(mcclust.ext)
source("R/dirichletprocess_custom_init.R")
source("R/density_clusterer.R")
source("R/ne_parts_pair_counting.R")
source("R/salso_custom.R")
source("R/credible_bounds.R")
source("R/rearrange_labels.R")
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
burn <- floor(nsims / 2)
for(d in 1:length(toy_datasets)){
  dataset_name <- names(toy_datasets)[d]
  if(!file.exists(paste0("output/toy_challenge/plot_obs_dpmm_", dataset_name, ".rds")) ||
     !file.exists(paste0("output/toy_challenge/plot_grid_dpmm_", dataset_name, ".rds"))){
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
    # Fit the DPMM  Model
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
    
    fn_samps_grid <- matrix(nrow = nsims, ncol = nrow(xy_grid))
    for(s in 1:nsims)
      fn_samps_grid[s,] <- PosteriorFunction(dp_mod, s)(xy_grid)
    plot_grid$f_pe <- colMeans(fn_samps_grid)
    rm(fn_samps_grid)
    
    # Get the Model-Based Cluster Allocations and WG Vertical Credible Bounds
    mixture_clustering_samps <- matrix(nrow = nsims, ncol = nrow(x))
    for(s in 1:nsims)
      mixture_clustering_samps[s,] <- dp_mod$labelsChain[[s]]
    # drop some burn-in samples
    mixture_clustering_samps <- mixture_clustering_samps[(burn + 1):nsims,]
    
    mixture_pe <- salso::salso(mixture_clustering_samps, 
                               maxZealousAttempts = 25)
    mixture_bounds <- credible_bounds(mixture_pe,
                                      mixture_clustering_samps)
    wg_mixture_bounds <- mcclust.ext::credibleball(mixture_pe, mixture_clustering_samps)
    plot_obs$mm_pe <- mixture_pe
    plot_obs$mm_vl <- mixture_bounds$min
    plot_obs$mm_vu <- mixture_bounds$max
    plot_obs$mm_vl_wg <- wg_mixture_bounds$c.lowervert[1,]
    plot_obs$mm_vu_wg <- wg_mixture_bounds$c.uppervert[1,]
    rm(dp_mod)
    
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
    
    saveRDS(plot_obs, paste0("output/toy_challenge/plot_obs_dpmm_", dataset_name, ".rds"))
    saveRDS(plot_grid, paste0("output/toy_challenge/plot_grid_dpmm_", dataset_name, ".rds"))
  }
}

# Load the Plotting Objects
plot_obs_tm <- readRDS("output/toy_challenge/plot_obs_dpmm_two_moons.rds")
plot_obs_nc <- readRDS("output/toy_challenge/plot_obs_dpmm_circles.rds")
plot_obs_ts <- readRDS("output/toy_challenge/plot_obs_dpmm_tsne.rds")

plot_grid_tm <- readRDS("output/toy_challenge/plot_grid_dpmm_two_moons.rds")
plot_grid_nc <- readRDS("output/toy_challenge/plot_grid_dpmm_circles.rds")
plot_grid_ts <- readRDS("output/toy_challenge/plot_grid_dpmm_tsne.rds")

# Rearrange labels to facilitate better color discrimination among the largest clusters
prep_labels <- function(plot_obs){
  label_columns <- setdiff(colnames(plot_obs), c('x', 'y'))
  for(lc in label_columns){
    plot_obs[[lc]] <- rearrange_labels(plot_obs[[lc]])
    plot_obs[[lc]] <- factor(plot_obs[[lc]], 
                           levels = 1:max(plot_obs[[lc]]))
  }
  return(plot_obs)
}
plot_obs_tm <- prep_labels(plot_obs_tm)
plot_obs_nc <- prep_labels(plot_obs_nc)
plot_obs_ts <- prep_labels(plot_obs_ts)

# Plot coloring palette creation function
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
color_vals_tm <- c(gg_color_hue(6), 
                   rep(RColorBrewer::brewer.pal(4, 'Greens'), 1000))
color_vals_ts <- c(gg_color_hue(10), 
                rep(RColorBrewer::brewer.pal(4, 'Greens'), 1000))

# Function to help us set plot titles based on cluster richness
ktitle <- function(x) {
  sprintf("K = %d (%d)", nlevels(x), sum(table(x) > 1))
}

# First we will visualize:
# - Density contours and heatplots
# - Model based cluster allocations
# - Density based cluster allocations
pf_tm <- ggplot(plot_grid_tm) +   
  geom_contour_filled(aes(x = x, y = y, z = f_pe)) + 
  guides(fill = 'none') + 
  labs(x = NULL, y = NULL)
pf_ts <- ggplot(plot_grid_ts) +   
  geom_contour_filled(aes(x = x, y = y, z = f_pe)) + 
  guides(fill = 'none') + 
  labs(x = NULL, y = NULL)

pmc_tm <- ggplot(plot_obs_tm) + 
  geom_point(aes(x = x, y = y, color = mm_pe), size = 0.5) + 
  stat_ellipse(aes(x = x, y = y, group = mm_pe)) + 
  scale_color_manual(values = color_vals_tm) +
  guides(color = 'none') +
  theme(plot.title = element_text(size=9)) +
  labs(title = ktitle(plot_obs_tm$mm_pe), x = NULL, y = NULL)
pmc_ts <- ggplot(plot_obs_ts) + 
  geom_point(aes(x = x, y = y, color = mm_pe), size = 0.25) + 
  stat_ellipse(aes(x = x, y = y, group = mm_pe)) + 
  scale_color_manual(values = color_vals_ts) +
  guides(color = 'none') + 
  theme(plot.title = element_text(size=9)) +
  labs(title = ktitle(plot_obs_ts$mm_pe), x = NULL, y = NULL)

pdc_tm <- ggplot() + 
  geom_point(aes(x = x, y = y, color = db_pe), 
             size = 0.5, data = plot_obs_tm) + 
  geom_contour(aes(x = x, y = y, z = f_pe), 
               breaks = 0.07, col = 'black',
               data = plot_grid_tm) + 
  scale_color_manual(values = color_vals_tm) +
  guides(color = 'none') + 
  theme(plot.title = element_text(size=9)) +
  labs(title = ktitle(plot_obs_tm$db_pe), x = NULL, y = NULL)
pdc_ts <- ggplot() + 
  geom_point(aes(x = x, y = y, color = db_pe), 
             size = 0.25, data = plot_obs_ts) + 
  ggforce::geom_mark_hull(
    aes(x = x, y = y, group = db_pe), 
    expand = 1e-2,
    radius = 1e-2,
    concavity = 2,
    data = plot_obs_ts %>% filter(db_pe != 0)) +
  scale_color_manual(values = color_vals_ts) +
  guides(color = 'none') + 
  theme(plot.title = element_text(size=9)) +
  labs(title = ktitle(plot_obs_ts$db_pe), x = NULL, y = NULL)

png("output/toy_challenge/compare_clusterings_dpmm.png", 
    width = 12, height = 8, units = 'in', res = 300)
grid.arrange(
  arrangeGrob(pf_tm, pf_ts, 
              top = textGrob(
                "DPMM Density Estimate",
                gp = gpar(fontface = 3, fontsize = 14)
              )),
  arrangeGrob(pmc_tm, pmc_ts, 
              top = textGrob(
                "Model-Based",
                gp = gpar(fontface = 3, fontsize = 14)
              )),
  arrangeGrob(pdc_tm, pdc_ts,
              top = textGrob(
                "BAND",
                gp = gpar(fontface = 3, fontsize = 14)
              )),
  top = textGrob("Clustering Point Estimates - DPMM", 
                 gp = gpar(fontsize = 18)),
  ncol = 3
)
dev.off()

# Then we will visualize model-based clusterings and 
# credible bounds for each dataset
pvl_tm <- ggplot(plot_obs_tm) + 
  geom_point(aes(x = x, y = y, color = mm_vl_wg), size = 0.5) + 
  stat_ellipse(aes(x = x, y = y, group = mm_vl_wg)) + 
  scale_color_manual(values = color_vals_tm) +
  guides(color = 'none') + 
  theme(plot.title = element_text(size=9)) +
  labs(title = ktitle(plot_obs_tm$mm_vl_wg), x = NULL, y = NULL)
pvl_ts <- ggplot(plot_obs_ts) + 
  geom_point(aes(x = x, y = y, color = mm_vl_wg), size = 0.25) + 
  stat_ellipse(aes(x = x, y = y, group = mm_vl_wg)) + 
  scale_color_manual(values = color_vals_ts) +
  guides(color = 'none') + 
  theme(plot.title = element_text(size=9)) +
  labs(title = ktitle(plot_obs_ts$mm_vl_wg), x = NULL, y = NULL)

pvu_tm <- ggplot(plot_obs_tm) + 
  geom_point(aes(x = x, y = y, color = mm_vu_wg), size = 0.5) + 
  stat_ellipse(aes(x = x, y = y, group = mm_vu_wg)) + 
  scale_color_manual(values = color_vals_tm) +
  guides(color = 'none') + 
  theme(plot.title = element_text(size=9)) +
  labs(title = ktitle(plot_obs_tm$mm_vu_wg), x = NULL, y = NULL)
pvu_ts <- ggplot(plot_obs_ts) + 
  geom_point(aes(x = x, y = y, color = mm_vu_wg), size = 0.25) + 
  stat_ellipse(aes(x = x, y = y, group = mm_vu_wg)) + 
  scale_color_manual(values = color_vals_ts) +
  guides(color = 'none') + 
  theme(plot.title = element_text(size=9)) +
  labs(title = ktitle(plot_obs_ts$mm_vu_wg), x = NULL, y = NULL)

png("output/toy_challenge/mbc_wg_bounds_dpmm.png", 
    width = 12, height = 8, units = 'in', res = 300)
grid.arrange(
  arrangeGrob(pvl_tm, pvl_ts, 
              top = textGrob(
                "Lower Bound",
                gp = gpar(fontface = 3, fontsize = 14)
              )),
  arrangeGrob(pmc_tm, pmc_ts, 
              top = textGrob(
                "Point Estimate",
                gp = gpar(fontface = 3, fontsize = 14)
              )),
  arrangeGrob(pvu_tm, pvu_ts,
              top = textGrob(
                "Upper Bound",
                gp = gpar(fontface = 3, fontsize = 14)
              )),
  top = textGrob("Model-Based Clustering - W&G Credible Balls", 
                 gp = gpar(fontsize = 18)),
  ncol = 3
)
dev.off()
