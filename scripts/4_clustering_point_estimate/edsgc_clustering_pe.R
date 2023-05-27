library(tidyverse)
source("R/ne_parts_pair_counting.R")


galaxy_data <- read.csv("data/intermediate_data/primary/edsgc.csv")
galaxy_clusters <- read.csv("data/intermediate_data/primary/edcci_clusters.csv")

clustering_samps <- readRDS("fitted_models/primary/edsgc_tree_dirichlet_clustering_samps.rds")

# We need to conserve memory and expedite optimization
# To do this, we first determine which observations are noise with very high probability
# In our optimization, we won't even consider alternate specifications for these, we will just leave them set at 0

always0_indcs <- which(apply(clustering_samps, 2, \(v) mean(v == 0) > 0.8))

clustering_samps <- clustering_samps[,-always0_indcs]

# Mimize Expected Loss
labels <- minimize_expected_loss(clustering_samps, 
                                 color_set = c(0,1))


labels_expanded <- rep(0, 41171)
labels_expanded[-always0_indcs] <- labels

labels_expanded <- factor(labels_expanded)
saveRDS(labels_expanded, file = "fitted_models/primary/edsgc_tree_dirichlet_clustering_pe.rds")

labels_expanded <- readRDS("fitted_models/primary/edsgc_tree_dirichlet_clustering_pe.rds")

labels <- as.numeric(as.character(labels_expanded))

label_counts <- table(labels)
substantial_clusters <- as.numeric(
  names(label_counts)[label_counts > 10])

masked_labels <- ifelse(labels %in% substantial_clusters, labels, 0)

labels <- factor(masked_labels)

labels_expanded <- labels

ggplot() + 
  geom_point(aes(x = x, y = y), color = 'grey', 
             alpha = 0.5, size = 0.1, data = galaxy_data[labels_expanded == 0,]) + 
  geom_point(aes(x = x, y = y, color = labels_expanded[labels_expanded != 0]), 
             alpha = 0.5, size = 0.1, data = galaxy_data[labels_expanded != 0,]) + 
  microViz::stat_chull(aes(x = x, y = y, color = labels_expanded[labels_expanded != 0]),
                       data = galaxy_data[labels_expanded != 0, ]) + 
  geom_point(aes(x = x, y = y), 
             size = 2, color = 'black', shape = 4, data = galaxy_clusters) + 
  guides(color = "none") + 
  # scale_x_continuous(expand = c(0,0)) + 
  # scale_y_continuous(expand = c(0,0)) + 
  labs(title = "EDSGC Data and EDCC Locations")

