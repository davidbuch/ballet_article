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
source("R/choose_lambda.R")
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

edcci <- as.matrix(edcci)
abell <- as.matrix(abell)

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

x.grid <- with(data, seq(min(x), max(x)*(1-1e-8), length.out = grid_res))
y.grid <- with(data, seq(min(y), max(y)*(1-1e-8), length.out = grid_res))

plot_grid <-   expand.grid(x = x.grid, y = y.grid)

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

z.grid <- matrix(NA, nrow=grid_res, ncol=grid_res)
for(i in seq_len(grid_res)) {
  for(j in seq_len(grid_res)) {
    z.grid[i,j] <- plot_grid$f_pe[i + (j-1)*grid_res]
  }
}


# The following be approximately equal to 
# with(data, 1/((max(x)-min(x))*(max(y)-min(y))))
avg_density <- mean(plot_grid$f_pe)
 
# The following code visualizes the density estimate and the level in 3D 
library(scales)

plot_persp <- function(levels, theta=60, phi=0, nlines=20, 
                       alpha.val=0.3, ...) {
  my_pal <- colorRampPalette(c("blue","red"))(length(levels))
  p <- persp(x.grid, y.grid, z.grid, theta=theta,
             phi=phi, ...)
  
  q.xs <- quantile(x.grid, seq_len(nlines)/nlines)
  
  for(i in seq_along(levels)) {
    for(l in seq_len(nlines)) {
      lines(trans3d(q.xs[l], 
                    y=y.grid, 
                    z=levels[i],  pmat=p), 
            col=alpha(my_pal[i], alpha.val)) 
    }
  }
}

levels <- c(1.8*avg_density, 
            2*avg_density,
            2.2*avg_density) 
{
  opar <- par()
  png("output/sky_survey_analysis/density_plot.png", 
    width = 5, height = 5, units = 'in', res = 1000)
  par(mfrow=c(2,3))
  for(i in 1:6) {
    plot_persp(levels, theta=60*i, 
               main=paste("Theta:", 60*i))
  }
  dev.off()
  par(opar)
}


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

# We will use scientific knowledge to compute the threshold lambda.
# Galaxies are known to form at overdensity = (f-average_f)/average_f >= 1.
# Thus this leads to the threshold of f >= 2 average_f
#must be around 1
overdensity_thresh <- seq(.8,1.2,by=0.1)

# The threshod index for which we will focus our analysis. Can be one of 1.. len(overdensity_thresh)
analysis_index <- length(overdensity_thresh)%/%2 + 1

# The lambda values corresponding to the overdensity threshold
lambda_scientific <- (1+ overdensity_thresh)*avg_density
target_quantiles <- ecdf(colMedians(f_samps))(lambda_scientific)
target_quantiles

BALLET_res <- rep(list(NULL), length(target_quantiles))

for(i in seq_along(target_quantiles)) {
  target_quantile <- target_quantiles[i]
  clustering_samps <- density_based_clusterer(X, f_samps, 
                                              cut_quantile = target_quantile)

  always0_indcs <- which(apply(clustering_samps, 2, \(v) mean(v == 0) > 0.8))
  non0_clustering_samps <- clustering_samps[,-always0_indcs]

  posterior_similarity_tensor <- compute_pst(non0_clustering_samps)
  posterior_difference_tensor <- compute_pdt(non0_clustering_samps)

  res <- list()
  res$pe <- salso_custom(clustering_samps,
                          pst = posterior_similarity_tensor,
                          pdt = posterior_difference_tensor,
                          always0_indcs = always0_indcs)

  bounds <- credible_ball_bounds_active_inactive(X, res$pe, clustering_samps)
  
  res$test_results <- data.frame(
      type = c(
               "BALLET_LB",
               "BALLET_PE",
               "BALLET_UB"
      ),
      sensitivity_edcci = c(
        sensitivity(X, res$lb, edcci),
        sensitivity(X, res$pe, edcci),
        sensitivity(X, res$ub, edcci)
      ),
      specificity_edcci = c(
        specificity(X, res$lb, edcci),
        specificity(X, res$pe, edcci),
        specificity(X, res$ub, edcci)
      ),
      exact_match_edcci = c(
        exact_match_frac(X, res$lb, edcci),
        exact_match_frac(X, res$pe, edcci),
        exact_match_frac(X, res$ub, edcci)
      ),
      sensitivity_abell = c(
        sensitivity(X, res$lb, abell),
        sensitivity(X, res$pe, abell),
        sensitivity(X, res$ub, abell)
      ),
      specificity_abell = c(
        specificity(X, res$lb, abell),
        specificity(X, res$pe, abell),
        specificity(X, res$ub, abell)
      ),
      exact_match_abell = c(
        exact_match_frac(X, res$lb, abell),
        exact_match_frac(X, res$pe, abell),
        exact_match_frac(X, res$ub, abell)
      )
  )
  
  BALLET_res[[i]] <- res
}

# Save the results to be able to reuse these.
names(BALLET_res) <- target_quantiles
saveRDS(BALLET_res, file = "BALLET_edsgc_res.rds")

# We will show these results for the median clustering.
rep_res <- BALLET_res[[analysis_index]]
target_quantile <- target_quantiles[[analysis_index]]
overdensity_choice <- overdensity_thresh[[analysis_index]]

data$pe <- rep_res$pe
data$ub <- rep_res$ub
data$lb <- rep_res$lb

# Next, we find the persistent clustering estimate.

#First form the cluster tree at various levels:
ctree <- do.call(cbind, map(BALLET_res, \(res) {
    #cf <- as.factor(res$pe)
    #nzlevels <- levels(cf) != "0"
    #levels(cf)[nzlevels] <- 1:sum(nzlevels)
    #cf
    res$pe
}))

colnames(ctree) <- paste("NoiseFrac", 
                          round(as.numeric(names(BALLET_res)), 
                                2))
clustree::clustree(ctree, prefix="NoiseFrac")
## TODO: save this plot.

per_clust <- select_persistent_clusters(ctree, prefix="NoiseFrac")
#per_clust <- as.factor(per_clust)
# assuming that 0 is the lowest level
#levels(per_clust) <- 0:(length(levels(per_clust))-1)

# Clean up the labels to that they take values 0, 1, \ldots, n.
# Note that label zero should always denote the noise label.
per_clust2 <- match(per_clust, as.integer(names(table(per_clust)))) - 1 
           
data$per <- per_clust

# Visualize the various BALLET results: 

# The alpha value to use for the datapoints.
bg_alpha <- 0.05

## First plot the lower bound.
enriched_data_lb <- enrich_small_clusters(data, 'lb', size_lb = 25)
p1 <- ggplot() + 
  geom_point(aes(x = x, y = y, color = factor(lb)), 
             size = 0.1, data = data %>% filter(lb != 0)) + 
  geom_point(aes(x = x, y = y), 
             size = 0.1, alpha = bg_alpha, color = 'black', 
             data = data %>% filter(lb == 0)) +
  stat_ellipse(aes(x = x, y = y, group = factor(lb)), 
               data = enriched_data_lb %>% filter(lb != 0)) +
  geom_point(aes(x = x, y = y, shape = catalogue), size = 2, data = clusters) +
  scale_shape_manual(name = "Catalogue",
                     labels = c("EDCCI", "Abell"),
                     values = c(4,3)) + 
  guides(color = 'none') + 
  labs(title = "BALLET 2.5%-ile Lower Bound",
       subtitle = paste("c=", overdensity_choice),
       x = NULL, y = NULL)

## Next, we plot the point estimate
enriched_data_pe <- enrich_small_clusters(data, 'pe', size_lb = 25)
p2 <- ggplot() + 
  geom_point(aes(x = x, y = y, color = factor(pe)),
             size = 0.1, data = data %>% filter(pe != 0)) + 
  geom_point(aes(x = x, y = y), size = 0.1, alpha = bg_alpha,
             color = 'black', data = data %>% filter(pe == 0)) +
  stat_ellipse(aes(x = x, y = y, group = factor(pe)), data = enriched_data_pe %>% filter(pe != 0)) +
  geom_point(aes(x = x, y = y, shape = catalogue), size = 2, data = clusters) +
  scale_shape_manual(name = "Catalogue",
                     labels = c("EDCCI", "Abell"),
                     values = c(4,3)) + 
  guides(color = 'none') + 
  labs(title = "BALLET Estimated Clusters",
       subtitle = paste("c=", overdensity_choice),
       x = NULL, y = NULL)

## Next, we plot the upper bound.
enriched_data_ub <- enrich_small_clusters(data, 'ub', size_lb = 25)
p3 <- ggplot() + 
  geom_point(aes(x = x, y = y, color = factor(ub)), 
             size = 0.1, data = data %>% filter(ub != 0)) + 
  geom_point(aes(x = x, y = y), size = 0.1, alpha = bg_alpha,
             color = 'black', data = data %>% filter(ub == 0)) +
  stat_ellipse(aes(x = x, y = y, group = factor(ub)), data = enriched_data_ub %>% filter(ub != 0)) +
  geom_point(aes(x = x, y = y, shape = catalogue), size = 2, data = clusters) +
  scale_shape_manual(name = "Catalogue",
                     labels = c("EDCCI", "Abell"),
                     values = c(4,3)) +
  guides(color = 'none') + 
  labs(title = "BALLET 97.5%-ile Upper Bound",
       subtitle = paste("c=", overdensity_choice),
       x = NULL, y = NULL)

# Finally, we plot the persistent clusterings
enriched_data_per <- enrich_small_clusters(data, 'per', size_lb = 25)
p4 <- ggplot() + 
  geom_point(aes(x = x, y = y, color = per), size = 0.1, data = data %>% filter(per != 0)) + 
  geom_point(aes(x = x, y = y), size = 0.1, alpha = bg_alpha, 
             color = 'black', data = data %>% filter(per == 0)) +
  stat_ellipse(aes(x = x, y = y, group = per), 
               data = enriched_data_per %>% filter(per != 0)) +
  geom_point(aes(x = x, y = y, shape = catalogue), size = 2, data = clusters) +
  scale_shape_manual(name = "Catalogue",
                     labels = c("EDCCI", "Abell"),
                     values = c(4,3)) + 
  guides(color = 'none') + 
  labs(title = "BALLET Persistent Clusters",
       subtitle = cat(c("c=", overdensity_thresh), sep=" "),
       x = NULL, y = NULL)


# Save the four plots to files
png("output/sky_survey_analysis/edsgc_ballet_pe.png",
    width = 5, height = 5, units = 'in', res = 300)
print(p2)
dev.off()

png("output/sky_survey_analysis/edsgc_ballet_bounds.png", 
    width = 10, height = 5, units = 'in', res = 300)
print(grid.arrange(p1, p3, nrow = 1))
dev.off()
  
png("output/sky_survey_analysis/edsgc_ballet_persistent.png",
    width = 5, height = 5, units = 'in', res = 300)
print(p4)
dev.off()

## Plot the valuations from older BALLET metrics
# map_dfr(BALLET_res, \(res) rownames_to_column(res$test_results),  
#         .id="noise_frac") %>%
# pivot_longer(
#   cols=starts_with("specificity")|starts_with("sensitivity")|starts_with("exact_match"),
#   names_to = "metric",
#   values_to = "score") %>%
#   mutate(catalogue=str_extract(metric, "(edcci|abell)"),
#          metric_type=str_extract(metric, "(specificity|sensitivity|exact_match)"),
#          noise_percent = ceiling(100*as.numeric(noise_frac))) -> df.ballet
#
# df.ballet %>%
#   filter(type == "BALLET_PE") %>%
# ggplot(aes(x=noise_percent, 
#              y=score, linetype=catalogue)) + 
#   geom_line() + 
#   facet_wrap(vars(metric_type), ncol=1) + 
#   ggtitle("BALLET p.e. vs existing catalogues.")
# ggsave("output/sky_survey_analysis/edsgc_ballet_vs_catalgoues.png")

# ----------------------------------------
# Find the DBSCAN clusters
# ----------------------------------------
# We will use the same minPts as we did for BALLET
target_quantile <- target_quantiles[[analysis_index]]
minPts <- ceiling(log2(nrow(X))) + 1
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
#minPts_NN_dists <- sort(minPts_NN_dists)
#eps <- minPts_NN_dists[round(napprox*(1 - target_quantile))]
epsilons <- quantile(minPts_NN_dists, probs = 1-target_quantiles)

DBSCAN_res <- map(epsilons, \(eps) {
 pe <- dbscan(X, eps = eps, minPts = minPts, borderPoints = FALSE)$cluster
 test_res <- data.frame(
      type = "DBSCAN",
      sensitivity_edcci = sensitivity(X, pe, edcci),
      specificity_edcci = specificity(X, pe, edcci),
      exact_match_edcci =  exact_match_frac(X, pe, edcci),
      sensitivity_abell =  sensitivity(X, pe, abell),
      specificity_abell =  specificity(X, pe, abell),
      exact_match_abell =  exact_match_frac(X, pe, abell)
      )
  list(clusters = pe, test_results=test_res)
  }
)



names(DBSCAN_res) <- target_quantiles

data$dbscan <- DBSCAN_res[[analysis_index]]$clusters

# map_dfr(DBSCAN_res, \(res) rownames_to_column(res$test_results),  
#         .id="noise_frac") %>%
#   pivot_longer(
#     cols=starts_with("specificity")|starts_with("sensitivity")|starts_with("exact_match"),
#     names_to = "metric",
#     values_to = "score") %>%
#   mutate(catalogue=str_extract(metric, "(edcci|abell)"),
#          metric_type=str_extract(metric, "(specificity|sensitivity|exact_match)"),
#          noise_percent = ceiling(100*as.numeric(noise_frac))) -> df.dbscan
#
# df.dbscan %>%
#   ggplot(aes(x=noise_percent, 
#              y=score, linetype=catalogue)) + 
#   geom_line() + 
#   facet_wrap(vars(metric_type), ncol=1) + 
#   ggtitle("DBSCAN  vs existing catalogues.")
#
#dbfit <-  dbscan(X, eps = eps, minPts = minPts, borderPoints = FALSE)$cluster
#data$dbscan <- dbfit$cluster


png("output/sky_survey_analysis/edsgc_dbscan_heuristic.png", 
    width = 5, height = 5, units = 'in', res = 300)
enriched_data <- enrich_small_clusters(data, 'dbscan', size_lb = 25)
ggplot() + 
  geom_point(aes(x = x, y = y, color = factor(dbscan)), 
             size = 0.1, data = data %>% filter(dbscan != 0)) + 
  geom_point(aes(x = x, y = y), 
             size = 0.1, alpha = bg_alpha, color = 'black', 
             data = data %>% filter(dbscan == 0)) +
  stat_ellipse(aes(x = x, y = y, group = factor(dbscan)), 
               data = enriched_data %>% filter(dbscan != 0)) +
  geom_point(aes(x = x, y = y, shape = catalogue), size = 2, data = clusters) +
  scale_shape_manual(name = "Catalogue",
                     labels = c("EDCCI", "Abell"),
                     values = c(4,3)) +
  guides(color = 'none') + 
  labs(title = "DBSCAN Estimated Clusters",
       subtitle = sprintf("MinPts: %d, Eps: %.2e", minPts, eps),
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
             size = 0.1, alpha = bg_alpha, color = 'black', 
             data = data %>% filter(dbscan_60 == 0)) +
  stat_ellipse(aes(x = x, y = y, group = factor(dbscan_60)), 
               data = enriched_data %>% filter(dbscan_60 != 0)) +
  geom_point(aes(x = x, y = y, shape = catalogue), size = 2, data = clusters) +
  scale_shape_manual(name = "Catalogue",
                     labels = c("EDCCI", "Abell"),
                     values = c(4,3)) +
  guides(color = 'none') + 
  labs(title = "DBSCAN Estimated Clusters",
       subtitle = sprintf("MinPts: %d, Eps: %.2e", minPts, eps),
       x = NULL, y = NULL)
dev.off()

## Now create a comparison table of all methods.
test_results <- data.frame(
  type = c("DBSCAN_60",
           "DBSCAN",
           "BALLET_LB",
           "BALLET_PE",
           "BALLET_UB",
           "BALLET (PER)"
  ),
  sensitivity_edcci = c(
    sensitivity(X, data$dbscan_60, edcci),
    sensitivity(X, data$dbscan,  edcci),
    sensitivity(X, data$lb, edcci),
    sensitivity(X, data$pe, edcci),
    sensitivity(X, data$ub, edcci),
    sensitivity(X, data$per, edcci)
  ),
  specificity_edcci = c(
    specificity(X, data$dbscan_60, edcci),
    specificity(X, data$dbscan, edcci),
    specificity(X, data$lb, edcci),
    specificity(X, data$pe, edcci),
    specificity(X, data$ub, edcci),
    specificity(X, data$per, edcci)
  ),
  exact_match_edcci = c(
    exact_match_frac(X, data$dbscan_60, edcci),
    exact_match_frac(X, data$dbscan, edcci),
    exact_match_frac(X, data$lb, edcci),
    exact_match_frac(X, data$pe, edcci),
    exact_match_frac(X, data$ub, edcci),
    exact_match_frac(X, data$per, edcci)
  ),
  sensitivity_abell = c(
    sensitivity(X, data$dbscan_60, abell),
    sensitivity(X, data$dbscan,  abell),
    sensitivity(X, data$lb, abell),
    sensitivity(X, data$pe, abell),
    sensitivity(X, data$ub, abell),
    sensitivity(X, data$per, abell)
  ),
  specificity_abell = c(
    specificity(X, data$dbscan_60, abell),
    specificity(X, data$dbscan, abell),
    specificity(X, data$lb, abell),
    specificity(X, data$pe, abell),
    specificity(X, data$ub, abell),
    specificity(X, data$per, abell)
  ),
  exact_match_abell = c(
    exact_match_frac(X, data$dbscan_60, abell),
    exact_match_frac(X, data$dbscan, abell),
    exact_match_frac(X, data$lb, abell),
    exact_match_frac(X, data$pe, abell),
    exact_match_frac(X, data$ub, abell),
    exact_match_frac(X, data$per, abell)
  )
)
saveRDS(test_results, 'output/sky_survey_analysis/edsgc_test_results.rds')
print(xtable(test_results), file = 'output/sky_survey_analysis/edsgc_test_results.txt')

