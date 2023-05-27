library(tidyverse)
set.seed(1234)
## Simulate Datasets

dir.create('data/clean_data', recursive = TRUE, showWarnings = FALSE)
dir.create('data/intermediate_data', recursive = TRUE, showWarnings = FALSE)

# The general idea for the collection of toy datasets
# comes from https://scikit-learn.org/stable/auto_examples/cluster/plot_cluster_comparison.html#sphx-glr-auto-examples-cluster-plot-cluster-comparison-py

# Generate datasets
nobs <- 500
spirals <- synthesis::data.gen.spirals(nobs, cycles = 1.5, s = 0.05, do.plot=FALSE)$x
noisy_circles <- synthesis::data.gen.circles(n = nobs, r_vec=c(1,1.5), start=runif(1,-1,1), s=0.05, do.plot=FALSE)$x
dstruct_moons <- fdm2id::data.twomoons(n = nobs, graph = FALSE, seed = 2)
two_moons <- cbind(dstruct_moons$X, dstruct_moons$Y)
blobs <- synthesis::data.gen.blobs(nobs = nobs, features=2, centers=5, sd=1, bbox=c(-10,10), do.plot=FALSE)$x
no_structure <- runif(nobs*2) %>% matrix(ncol = 2)

# anisotropic data
set.seed(3)
pre_X <- synthesis::data.gen.blobs(nobs = nobs, features=2, centers=3, sd=0.5, bbox=c(-6,6), do.plot=FALSE)$x
transformation <- matrix(c(0.6,-0.6,
                           -0.4, 0.8), 
                         nrow = 2, byrow = TRUE)
aniso <- pre_X %*% transformation

# blobs with varied variances
set.seed(0)
dstruct <- synthesis::data.gen.blobs(nobs = nobs, features=2, centers=3, sd=1, bbox=c(-6,6), do.plot=FALSE)
varied <- dstruct$x; bclasses <- dstruct$classes 
bclass_means <- t(sapply(bclasses, \(cl) colMeans(varied[bclasses == cl,])))
bclass_sds <- sapply(bclasses, \(cl) c(1.0, 0.5, 2.5)[cl])
varied <- (varied - bclass_means)*bclass_sds + bclass_means


simulated_datasets <- list(spirals = spirals,
                           two_moons = two_moons, 
                           blobs = blobs, 
                           aniso = aniso, 
                           varied = varied, 
                           no_structure = no_structure,
                           noisy_circles = noisy_circles)

simulated_datasets <- lapply(simulated_datasets, 
                             \(mat) data.frame(x = scale(mat[,1]), 
                                               y = scale(mat[,2])))

for(d in 1:length(simulated_datasets)){
  dataset_name <- names(simulated_datasets)[d]
  data <- simulated_datasets[[d]]
  
  write.csv(data, file = paste0("data/intermediate_data/", dataset_name,".csv"), row.names = FALSE)
  saveRDS(data, file = paste0("data/clean_data/", dataset_name,".rds"))
}

# tSNE Data
# Toy data set 'tSNE' comes from
# https://www.reneshbedre.com/blog/dbscan-python.html
# if we use this in the article we should reach out to him
# for permission
tsne <- read.csv("https://reneshbedre.github.io/assets/posts/tsne/tsne_scores.csv", 
                 check.names = FALSE)
tsne <- tsne %>% 
  mutate(across(where(is.numeric), ~ as.numeric(scale(.)))) %>%
  rename(x = `t-SNE-1`, y = `t-SNE-2`)
write.csv(tsne, file = "data/intermediate_data/tsne.csv", row.names = FALSE)
saveRDS(tsne, file = "data/clean_data/tsne.rds")



# library(palmerpenguins)

# # 1. Roeder Galaxy data
# write(MASS::galaxies, file = "data/intermediate_data/toy/galaxies.csv")
# saveRDS(MASS::galaxies, file = "data/clean_data/toy/galaxies.rds")

# # 2. Palmer Penguin data
# penguins <- penguins %>% 
#   dplyr::select(bill_length_mm, 
#          bill_depth_mm, 
#          flipper_length_mm, 
#          body_mass_g,
#          species) %>%
#   drop_na()
# penguins_labels <- penguins %>% 
#   pull(species)
# penguins_data <- penguins %>% 
#   dplyr::select(-species)
# write.csv(penguins_data, file = "data/intermediate_data/toy/penguins.csv", row.names = FALSE)
# saveRDS(penguins_data, "data/clean_data/toy/penguins.rds")
# saveRDS(penguins_labels, "data/clean_data_labels/toy/penguins.rds")


# # 3. UCI Wine data
# wine <- read.csv("data/raw_data/wine.data", header = FALSE, 
#                  col.names = c("Class", "Alcohol", "Malic Acid", 
#                                "Ash", "Alcalinity", "Magnesium", 
#                                "Phenols", "Flavanoids", "Nonflavanoids", 
#                                "Proanthocyanins", "Color", "Hue", 
#                                "OD_Ratio", "Proline")
#                  )
# wine <- wine %>% mutate(Class = factor(Class))
# wine_labels <- wine %>%
#   pull(Class)
# wine_data <- wine %>%
#   select(-Class)
# write.csv(wine_data, file = "data/intermediate_data/toy/wine.csv", row.names = FALSE)
# saveRDS(wine_data, "data/clean_data/toy/wine.rds")
# saveRDS(wine_labels, "data/clean_data_labels/toy/wine.rds")




