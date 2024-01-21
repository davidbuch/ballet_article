library(tidyverse)
library(gridExtra)
library(dbscan)
library(xtable)
source("R/random_histogram_model.R")
source("R/density_clusterer.R")
source("R/ne_parts_pair_counting.R")
source("R/salso_custom.R")
source("R/credible_bounds.R")
source("R/simstudy_accuracy_measures.R")
dir.create("output/sky_survey_analysis", showWarnings = FALSE, recursive = TRUE)

random_seed <- 1234
set.seed(random_seed)

# ----------------------------------------
# Load the Dataset
# ----------------------------------------
data <- read.csv("data/intermediate_data/edsgc.csv")
X <- as.matrix(data)
edcci <- read.csv("data/intermediate_data/edcci_clusters.csv")
abell <- read.csv("data/intermediate_data/abell_clusters.csv")

clusters <- rbind(data.frame(edcci, catalogue = "EDCCI"),
                  data.frame(abell, catalogue = "Abell"))

png("output/sky_survey_analysis/edsgc_data.png", 
    width = 5, height = 5, units = 'in', res = 300)
ggplot() +
  geom_point(aes(x = x, y = y), alpha = 0.2, size = 0.1, color = "grey50", data = data) +
  geom_point(aes(x = x, y = y, shape = catalogue), size = 2, data = clusters) +
  scale_shape_manual(name = "Cluster\nCatalogue",
                     labels = c("EDCCI", "Abell"),
                     values = c(4,3)) +
  guides(fill = 'none') + 
  labs(title = "EDSGC Galaxies and Catalogue Locations", 
       x = NULL, y = NULL)
dev.off()

# ----------------------------------------
# Fit the Random Histogram Model
# ----------------------------------------
nobs <- nrow(data)
nsamps <- 1000

prior_alpha <- 1
res <- 50
nhists <- 50

hist_data <- get_hist_data(data, res, nhists)

# Plot the Density Estimate on a Grid
grid_res <- 100
plot_grid <- with(data, 
                  expand.grid(x = seq(min(x), max(x)*(1-1e-8), length.out = grid_res),
                              y = seq(min(y), max(y)*(1-1e-8), length.out = grid_res))
)
plot_grid$f_pe <- rep(0, nrow(plot_grid))
pb <- txtProgressBar(width = 30, style = 3)
for(s in 1:nsamps){
  tryCatch({
    plot_grid$f_pe <- plot_grid$f_pe + 
      sample_density(plot_grid,
                     prior_alpha,
                     hist_data
      )
  }, error = function(e){s <- s - 1})
  setTxtProgressBar(pb, s/nsamps)
}
close(pb)
plot_grid$f_pe <- plot_grid$f_pe / nsamps


png("output/sky_survey_analysis/edsgc_randhist_density.png", 
    width = 5, height = 5, units = 'in', res = 300)
ggplot() + 
  geom_contour_filled(aes(x = x, y = y, z = log(f_pe)), data = plot_grid) + 
  geom_point(aes(x = x, y = y, shape = catalogue), size = 2, data = clusters) +
  scale_shape_manual(name = "Cluster\nCatalogue",
                     labels = c("EDCCI", "Abell"),
                     values = c(4,3)) +
  guides(fill = 'none') + 
  labs(title = "Log of Posterior Expected Density", 
                               x = NULL, y = NULL)
dev.off()
rm(plot_grid)

# Sample from the posterior of the density function
f_samps <- matrix(nrow = nsamps, ncol = nrow(data))
pb <- txtProgressBar(width = 30, style = 3)
for(s in 1:nsamps){
  tryCatch({
    f_samps[s,] <- sample_density(data,
                                  prior_alpha,
                                  hist_data)
  }, error = function(e){s <- s - 1})
  setTxtProgressBar(pb, s/nsamps)
}
close(pb)



# ----------------------------------------
# Find the BALLET clusters
# ----------------------------------------
# We will use a scientifically motivated lambda rather than a target_quantile

# target_quantile <- 0.9
# 
# approximate_number_of_clusters <- 100
# points_per_cluster <- (1 - target_quantile)*nobs / 
#   approximate_number_of_clusters
# minPts <- round(points_per_cluster / 4)
# print(sprintf("Hueristic Suggested minPts: %d", minPts))
# # Stick with the default based on simulation analysis
# split_err_prob <- 0.05
# 
# clustering_samps <- density_based_clusterer(X, 
#                                             f_samps, 
#                                             cut_quantile = target_quantile,
#                                             minPts = minPts,
#                                             split_err_prob = split_err_prob)


# ----------------------------------------
# Find the BALLET clusters (SPEP 0.01)
# ----------------------------------------
# Scientifically motivated lambda: overdensity exceeds 1
# - equivalent to 2 times the average density on the space
# - equivalent to 2 times (1 / Area)
lambda <- 2 / prod(apply(plot_grid[c('x','y')], 2, max))
delta <- compute_delta(data, f_samps, lambda, 10, 0.01)

clustering_samps <- density_based_clusterer(X,
                                            f_samps,
                                            lambda_delta = c(lambda, delta))

always0_indcs <- which(apply(clustering_samps, 2, \(v) mean(v == 0) > 0.8))
non0_clustering_samps <- clustering_samps[,-always0_indcs]

posterior_similarity_tensor <- compute_pst(non0_clustering_samps)
posterior_difference_tensor <- compute_pdt(non0_clustering_samps)

data$pe <- salso_custom(clustering_samps,
                        pst = posterior_similarity_tensor,
                        pdt = posterior_difference_tensor,
                        always0_indcs = always0_indcs)

bounds <- credible_ball_bounds_active_inactive(X, data$pe, clustering_samps)
data$lb <- bounds$lower
data$ub <- bounds$upper

enriched_data_lb <- enrich_small_clusters(data, 'lb', size_lb = 25)
p1 <- ggplot() + 
  geom_point(aes(x = x, y = y, color = factor(lb)), 
             size = 0.1, data = data %>% filter(lb != 0)) + 
  geom_point(aes(x = x, y = y), 
             size = 0.1, alpha = 0.2, color = 'black', 
             data = data %>% filter(lb == 0)) +
  stat_ellipse(aes(x = x, y = y, group = factor(lb)), 
               data = enriched_data_lb %>% filter(lb != 0)) +
  geom_point(aes(x = x, y = y, shape = catalogue), size = 2, data = clusters) +
  scale_shape_manual(name = "Cluster\nCatalogue",
                     labels = c("EDCCI", "Abell"),
                     values = c(4,3)) + 
  guides(color = 'none') + 
  labs(title = "BALLET 2.5%-ile Lower Bound",
       subtitle = sprintf("minPts: %d, spep: %.2f", minPts, split_err_prob),
       x = NULL, y = NULL)

enriched_data_pe <- enrich_small_clusters(data, 'pe', size_lb = 25)
p2 <- ggplot() + 
  geom_point(aes(x = x, y = y, color = factor(pe)), size = 0.1, data = data %>% filter(pe != 0)) + 
  geom_point(aes(x = x, y = y), size = 0.1, alpha = 0.2, color = 'black', data = data %>% filter(pe == 0)) +
  stat_ellipse(aes(x = x, y = y, group = factor(pe)), data = enriched_data_pe %>% filter(pe != 0)) +
  geom_point(aes(x = x, y = y, shape = catalogue), size = 2, data = clusters) +
  scale_shape_manual(name = "Cluster\nCatalogue",
                     labels = c("EDCCI", "Abell"),
                     values = c(4,3)) + 
  guides(color = 'none') + 
  labs(title = "BALLET Estimated Clusters and Catalogue Locations",
       subtitle = sprintf("minPts: %d, spep: %.2f", minPts, split_err_prob),
       x = NULL, y = NULL)

enriched_data_ub <- enrich_small_clusters(data, 'ub', size_lb = 25)
p3 <- ggplot() + 
  geom_point(aes(x = x, y = y, color = factor(ub)), size = 0.1, data = data %>% filter(ub != 0)) + 
  geom_point(aes(x = x, y = y), size = 0.1, alpha = 0.2, color = 'black', data = data %>% filter(ub == 0)) +
  stat_ellipse(aes(x = x, y = y, group = factor(ub)), data = enriched_data_ub %>% filter(ub != 0)) +
  geom_point(aes(x = x, y = y, shape = catalogue), size = 2, data = clusters) +
  scale_shape_manual(name = "Cluster\nCatalogue",
                     labels = c("EDCCI", "Abell"),
                     values = c(4,3)) +
  guides(color = 'none') + 
  labs(title = "BALLET 97.5%-ile Upper Bound",
       subtitle = sprintf("minPts: %d, spep: %.2f", minPts, split_err_prob),
       x = NULL, y = NULL)

png("output/sky_survey_analysis/edsgc_ballet_pe_01.png", 
    width = 5, height = 5, units = 'in', res = 300)
print(p2)
dev.off()

png("output/sky_survey_analysis/edsgc_ballet_bounds_01.png", 
    width = 10, height = 5, units = 'in', res = 300)
print(grid.arrange(p1, p3, nrow = 1))
dev.off()

# ----------------------------------------
# Find the BALLET clusters (SPEP 0.05)
# ----------------------------------------
# Scientifically motivated lambda: overdensity exceeds 1
# - equivalent to 2 times the average density on the space
# - equivalent to 2 times (1 / Area)
lambda <- 2 / prod(apply(plot_grid[c('x','y')], 2, max))
delta <- compute_delta(data, f_samps, lambda, 10, 0.05)

clustering_samps <- density_based_clusterer(X,
                                            f_samps,
                                            lambda_delta = c(lambda, delta))

always0_indcs <- which(apply(clustering_samps, 2, \(v) mean(v == 0) > 0.8))
non0_clustering_samps <- clustering_samps[,-always0_indcs]

posterior_similarity_tensor <- compute_pst(non0_clustering_samps)
posterior_difference_tensor <- compute_pdt(non0_clustering_samps)

data$pe <- salso_custom(clustering_samps,
                        pst = posterior_similarity_tensor,
                        pdt = posterior_difference_tensor,
                        always0_indcs = always0_indcs)

bounds <- credible_ball_bounds_active_inactive(X, data$pe, clustering_samps)
data$lb <- bounds$lower
data$ub <- bounds$upper

enriched_data_lb <- enrich_small_clusters(data, 'lb', size_lb = 25)
p1 <- ggplot() + 
  geom_point(aes(x = x, y = y, color = factor(lb)), 
             size = 0.1, data = data %>% filter(lb != 0)) + 
  geom_point(aes(x = x, y = y), 
             size = 0.1, alpha = 0.2, color = 'black', 
             data = data %>% filter(lb == 0)) +
  stat_ellipse(aes(x = x, y = y, group = factor(lb)), 
               data = enriched_data_lb %>% filter(lb != 0)) +
  geom_point(aes(x = x, y = y, shape = catalogue), size = 2, data = clusters) +
  scale_shape_manual(name = "Cluster\nCatalogue",
                     labels = c("EDCCI", "Abell"),
                     values = c(4,3)) + 
  guides(color = 'none') + 
  labs(title = "BALLET 2.5%-ile Lower Bound",
       subtitle = sprintf("minPts: %d, spep: %.2f", minPts, split_err_prob),
       x = NULL, y = NULL)

enriched_data_pe <- enrich_small_clusters(data, 'pe', size_lb = 25)
p2 <- ggplot() + 
  geom_point(aes(x = x, y = y, color = factor(pe)), size = 0.1, data = data %>% filter(pe != 0)) + 
  geom_point(aes(x = x, y = y), size = 0.1, alpha = 0.2, color = 'black', data = data %>% filter(pe == 0)) +
  stat_ellipse(aes(x = x, y = y, group = factor(pe)), data = enriched_data_pe %>% filter(pe != 0)) +
  geom_point(aes(x = x, y = y, shape = catalogue), size = 2, data = clusters) +
  scale_shape_manual(name = "Cluster\nCatalogue",
                     labels = c("EDCCI", "Abell"),
                     values = c(4,3)) + 
  guides(color = 'none') + 
  labs(title = "BALLET Estimated Clusters and Catalogue Locations",
       subtitle = sprintf("minPts: %d, spep: %.2f", minPts, split_err_prob),
       x = NULL, y = NULL)

enriched_data_ub <- enrich_small_clusters(data, 'ub', size_lb = 25)
p3 <- ggplot() + 
  geom_point(aes(x = x, y = y, color = factor(ub)), size = 0.1, data = data %>% filter(ub != 0)) + 
  geom_point(aes(x = x, y = y), size = 0.1, alpha = 0.2, color = 'black', data = data %>% filter(ub == 0)) +
  stat_ellipse(aes(x = x, y = y, group = factor(ub)), data = enriched_data_ub %>% filter(ub != 0)) +
  geom_point(aes(x = x, y = y, shape = catalogue), size = 2, data = clusters) +
  scale_shape_manual(name = "Cluster\nCatalogue",
                     labels = c("EDCCI", "Abell"),
                     values = c(4,3)) +
  guides(color = 'none') + 
  labs(title = "BALLET 97.5%-ile Upper Bound",
       subtitle = sprintf("minPts: %d, spep: %.2f", minPts, split_err_prob),
       x = NULL, y = NULL)

png("output/sky_survey_analysis/edsgc_ballet_pe_05.png", 
    width = 5, height = 5, units = 'in', res = 300)
print(p2)
dev.off()

png("output/sky_survey_analysis/edsgc_ballet_bounds_05.png", 
    width = 10, height = 5, units = 'in', res = 300)
print(grid.arrange(p1, p3, nrow = 1))
dev.off()




# ----------------------------------------
# Find the DBSCAN clusters
# ----------------------------------------
# We will use the same minPts as we did for BALLET
tX <- t(X)
napprox <- 1000 # find epsilon based on an approximation of the full data
minPts_NN_dists <- c()
row_samps <- sample(1:nobs, napprox)
for(ridx in row_samps){
  x <- tX[,ridx]
  dists <- sqrt(colSums((tX - x)^2))
  minPts_NN_dists <- c(minPts_NN_dists,
                       sort(dists, partial = minPts)[minPts])
}
minPts_NN_dists <- sort(minPts_NN_dists)
eps <- minPts_NN_dists[round(napprox*(1 - target_quantile))]

dbfit <- dbscan(X, eps = eps, minPts = minPts, borderPoints = FALSE)
data$dbscan <- dbfit$cluster


png("output/sky_survey_analysis/edsgc_dbscan_heuristic.png", 
    width = 5, height = 5, units = 'in', res = 300)
enriched_data <- enrich_small_clusters(data, 'dbscan', size_lb = 25)
ggplot() + 
  geom_point(aes(x = x, y = y, color = factor(dbscan)), 
             size = 0.1, data = data %>% filter(dbscan != 0)) + 
  geom_point(aes(x = x, y = y), 
             size = 0.1, alpha = 0.2, color = 'black', 
             data = data %>% filter(dbscan == 0)) +
  stat_ellipse(aes(x = x, y = y, group = factor(dbscan)), 
               data = enriched_data %>% filter(dbscan != 0)) +
  geom_point(aes(x = x, y = y, shape = catalogue), size = 2, data = clusters) +
  scale_shape_manual(name = "Cluster\nCatalogue",
                     labels = c("EDCCI", "Abell"),
                     values = c(4,3)) +
  guides(color = 'none') + 
  labs(title = "DBSCAN Estimated Clusters and Catalogue Locations",
       subtitle = sprintf("minPts: %d, eps: %.2e", minPts, eps),
       x = NULL, y = NULL)
dev.off()



# ----------------------------------------
# Find the DBSCAN clusters II - minPts = 60
# ----------------------------------------
# We now try minPts as was calibrated during simulation
minPts <- 60

napprox <- 1000 # find epsilon based on an approximation of the full data
minPts_NN_dists <- c()
row_samps <- sample(1:nobs, napprox)
for(ridx in row_samps){
  x <- tX[,ridx]
  dists <- sqrt(colSums((tX - x)^2))
  minPts_NN_dists <- c(minPts_NN_dists,
                       sort(dists, partial = minPts)[minPts])
}
minPts_NN_dists <- sort(minPts_NN_dists)
eps <- minPts_NN_dists[round(napprox*(1 - target_quantile))]

dbfit <- dbscan(X, eps = eps, minPts = minPts, borderPoints = FALSE)
data$dbscan_60 <- dbfit$cluster


png("output/sky_survey_analysis/edsgc_dbscan_60.png", 
    width = 5, height = 5, units = 'in', res = 300)
enriched_data <- enrich_small_clusters(data, 'dbscan_60', size_lb = 25)
ggplot() + 
  geom_point(aes(x = x, y = y, color = factor(dbscan_60)), 
             size = 0.1, data = data %>% filter(dbscan_60 != 0)) + 
  geom_point(aes(x = x, y = y), 
             size = 0.1, alpha = 0.2, color = 'black', 
             data = data %>% filter(dbscan_60 == 0)) +
  stat_ellipse(aes(x = x, y = y, group = factor(dbscan_60)), 
               data = enriched_data %>% filter(dbscan_60 != 0)) +
  geom_point(aes(x = x, y = y, shape = catalogue), size = 2, data = clusters) +
  scale_shape_manual(name = "Cluster\nCatalogue",
                     labels = c("EDCCI", "Abell"),
                     values = c(4,3)) +
  guides(color = 'none') + 
  labs(title = "DBSCAN Estimated Clusters and Catalogue Locations",
       subtitle = sprintf("minPts: %d, eps: %.2e", minPts, eps),
       x = NULL, y = NULL)
dev.off()

edcci <- as.matrix(edcci)
abell <- as.matrix(abell)

test_results <- data.frame(
  type = c("DBSCAN_60",
           "DBSCAN",
           "BALLET_LB",
           "BALLET_PE",
           "BALLET_UB"
  ),
  sensitivity_edcci = c(
    sensitivity(X, data$dbscan_60, edcci),
    sensitivity(X, data$dbscan,  edcci),
    sensitivity(X, data$lb, edcci),
    sensitivity(X, data$pe, edcci),
    sensitivity(X, data$ub, edcci)
  ),
  specificity_edcci = c(
    specificity(X, data$dbscan_60, edcci),
    specificity(X, data$dbscan, edcci),
    specificity(X, data$lb, edcci),
    specificity(X, data$pe, edcci),
    specificity(X, data$ub, edcci)
  ),
  sensitivity_abell = c(
    sensitivity(X, data$dbscan_60, abell),
    sensitivity(X, data$dbscan,  abell),
    sensitivity(X, data$lb, abell),
    sensitivity(X, data$pe, abell),
    sensitivity(X, data$ub, abell)
  ),
  specificity_abell = c(
    specificity(X, data$dbscan_60, abell),
    specificity(X, data$dbscan, abell),
    specificity(X, data$lb, abell),
    specificity(X, data$pe, abell),
    specificity(X, data$ub, abell)
  )
)
saveRDS(test_results, 'output/sky_survey_analysis/edsgc_test_results.rds')
print(xtable(test_results), file = 'output/sky_survey_analysis/edsgc_test_results.txt')

