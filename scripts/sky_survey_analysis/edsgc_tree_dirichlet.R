library(tidyverse)

galaxy_data <- read.csv("data/intermediate_data/primary/edsgc.csv")
galaxy_clusters <- read.csv("data/intermediate_data/primary/edcci_clusters.csv")

probs <- fit_tree_dirichlet(galaxy_data)
saveRDS(probs, file = "fitted_models/primary/edsgc_tree_dirichlet.rds")
