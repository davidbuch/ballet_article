library(tidyverse)
library(grid)
library(gridExtra)
library(dirichletprocess)
source("R/density_clusterer.R")
source("R/salso_custom.R")
source("R/ne_parts_pair_counting.R")
source("R/rearrange_labels.R")

# Unrelated to the eventual def. of quantile bounds, but choice can improve or 
# deteriorate the clarity of the bounds' properties.
DENSITY_THRESHOLD_QUANTILE <- 0.125
SPLIT_ERR_PROB <- 0.01

# dp_mod <- readRDS("../bad_clustering_article/fitted_models/toy_data/tsne_dpmm_c1.rds")
# x <- dp_mod$data
# nsims <- length(dp_mod$alphaChain)

# Load Data and Set Parameters
x <- scale(as.matrix(readRDS("data/clean_data/tsne.rds")))
nobs <- nrow(x)
nsims <- 1000
alpha <- 0.05

# Fit the DPMM  Model
g0Priors <- list(mu0 = rep(0,ncol(x)),
                 Lambda = diag(ncol(x)),
                 kappa0 = 0.1, #0.01,
                 nu = 10) # 200
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

# NOTE: WE'RE TOYING WITH THE DENSITY THRESHOLD QUANTILE
density_clustering_samps <- 
  density_based_clusterer(x, fn_samps_obs, 
                          cut_quantile = DENSITY_THRESHOLD_QUANTILE,
                          split_err_prob = SPLIT_ERR_PROB)
density_pe <- salso_custom(density_clustering_samps)

# Method 1 - noise/non-noise level bounds
source('R/credible_bounds.R')
bounds <- credible_bounds_active_inactive(x, density_clustering_samps)
active_lb <- bounds$lower
active_ub <- bounds$upper

# Method 2 - whatever we were doing before
bounds <- credible_bounds(density_pe, density_clustering_samps)
cblb <- bounds$min
cbub <- bounds$max

# Other computations of Interest
S <- nrow(density_clustering_samps)
clustering_distances <- vector(length = S)
pb <- txtProgressBar(width = 30, style = 3)
for(s in 1:S){
  clustering_distances[s] <- ne_loss(density_pe, 
                                         density_clustering_samps[s,])
  setTxtProgressBar(pb,value = s/S)
}
close(pb)

epsilon_bound <- quantile(clustering_distances, 0.95, names = FALSE)

ne_loss(density_pe, active_ub)/epsilon_bound
ne_loss(density_pe, active_lb)/epsilon_bound
ne_loss(density_pe, cbub)/epsilon_bound
ne_loss(density_pe, cblb)/epsilon_bound


# Helper functions
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
nunique <- function(lbs){
  length(unique(lbs))
}
color_vals_ts <- c(gg_color_hue(10), 
                   rep(RColorBrewer::brewer.pal(4, 'Greens'), 1000))

density_pe <- rearrange_labels(density_pe)
density_pe <- replace(density_pe, density_pe == 0, NA)

active_lb <- rearrange_labels(active_lb)
active_ub <- rearrange_labels(active_ub)

cblb <- rearrange_labels(cblb)
cbub <- rearrange_labels(cbub)

active_lb <- replace(active_lb, active_lb == 0, NA)
active_ub <- replace(active_ub, active_ub == 0, NA)

cblb <- replace(cblb, cblb == 0, NA)
cbub <- replace(cbub, cbub == 0, NA)

p1 <- ggplot() + 
  geom_point(aes(x = x[,1], y = x[,2], color = factor(density_pe)), size = 0.1) + 
  guides(color = 'none') + scale_color_manual(values = color_vals_ts) +
  labs(title = sprintf("Point Estimate (%d)", nunique(density_pe)), 
       x = NULL, y = NULL)
p2 <- ggplot() + 
  geom_point(aes(x = x[,1], y = x[,2], color = factor(active_lb)), size = 0.1) + 
  guides(color = 'none') + scale_color_manual(values = color_vals_ts) +
  labs(title = sprintf("Active Lower Bound (%d)", nunique(active_lb)), 
       x = NULL, y = NULL)
p3 <- ggplot() + 
  geom_point(aes(x = x[,1], y = x[,2], color = factor(active_ub)), size = 0.1) + 
  guides(color = 'none') + scale_color_manual(values = color_vals_ts) + 
  labs(title = sprintf("Active Upper Bound (%d)", nunique(active_ub)), 
       x = NULL, y = NULL)

p4 <- ggplot() + 
  geom_point(aes(x = x[,1], y = x[,2], color = factor(cblb)), size = 0.1) + 
  guides(color = 'none') + scale_color_manual(values = color_vals_ts) +
  labs(title = sprintf("Lattice Min Lower Bound (%d)", nunique(cblb)), 
       x = NULL, y = NULL)
p5 <- ggplot() + 
  geom_point(aes(x = x[,1], y = x[,2], color = factor(cbub)), size = 0.1) + 
  guides(color = 'none') + scale_color_manual(values = color_vals_ts) + 
  labs(title = sprintf("Lattice Max Upper Bound (%d)", nunique(cbub)), 
       x = NULL, y = NULL)

png("output/toy_challenge/density_bound_analysis.png", 
    width = 15, height = 10, units = 'in', res = 300)
grid.arrange(p2, p1, p3,
             p4, p1, p5,
             nrow = 2, ncol = 3, 
             top = textGrob("Bound Analysis - Density Based Clustering", 
                            gp = gpar(fontsize = 18)))
dev.off()

