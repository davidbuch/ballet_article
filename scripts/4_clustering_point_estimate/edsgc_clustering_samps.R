library(tidyverse)
source("R/density_clusterer.R")
source("R/ne_parts_pair_counting.R")

galaxy_data <- read.csv("data/intermediate_data/primary/edsgc.csv")
galaxy_clusters <- read.csv("data/intermediate_data/primary/edcci_clusters.csv")
fsamps <- readRDS("fitted_models/primary/edsgc_tree_dirichlet.rds")

# Clustering 
lambda <- 0.9
clustering_samps <- density_based_clusterer(galaxy_data, 
                                  fsamps, 
                                  cut_quantile = lambda,
                                  minPts = 2,
                                  split_err_prob = 0.10
)
saveRDS(clustering_samps, "fitted_models/primary/edsgc_tree_dirichlet_clustering_samps.rds")

rm(galaxy_data, fsamps)

# Mimize Expected Loss
labels <- minimize_expected_loss(clustering_samps)

labels <- factor(clustering_samps[1000,])

ggplot() + 
  geom_point(aes(x = x, y = y, color = probs > quantile(probs, lambda)), 
             alpha = 0.5, size = 0.1, data = galaxy_data) + 
  geom_point(aes(x = x, y = y), 
             size = 0.1, color = "black", data = galaxy_clusters) + 
  guides(color = "none") + 
  # scale_x_continuous(expand = c(0,0)) + 
  # scale_y_continuous(expand = c(0,0)) + 
  labs(title = "EDSGC Data and EDCC Locations")

ggplot() + 
  geom_point(aes(x = x, y = y), color = 'grey', 
             alpha = 0.5, size = 0.1, data = galaxy_data[labels == 0,]) + 
  geom_point(aes(x = x, y = y, color = labels[labels != 0]), 
             alpha = 0.5, size = 0.1, data = galaxy_data[labels != 0,]) + 
  microViz::stat_chull(aes(x = x, y = y, color = labels[labels != 0]),
                       data = galaxy_data[labels != 0, ]) + 
  geom_point(aes(x = x, y = y), 
             size = 0.1, color = 'black', data = galaxy_clusters) + 
  guides(color = "none") + 
  # scale_x_continuous(expand = c(0,0)) + 
  # scale_y_continuous(expand = c(0,0)) + 
  labs(title = "EDSGC Data and EDCC Locations")
