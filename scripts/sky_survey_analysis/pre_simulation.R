library(tidyverse)
library(gridExtra)
library(dbscan)
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
# Simulate the Dataset
# ----------------------------------------
N <- 40000 # n obs
R <- 0.90 # noise fraction
K <- 42 # number of non-noise clusters

# simulate active allocations from a dirichlet dsn
active_parts <- rgamma(K, 3)
active_parts <- active_parts / sum(active_parts)

png("output/sky_survey_analysis/simulation_component_weights.png", 
    width = 5, height = 5, units = 'in', res = 300)
barplot((1 - R)*active_parts, 
        xlab = 'Component ID', 
        ylab = 'Assignment Probabilities', 
        main = "Non-Noise Component Weights")
dev.off()

pi <- c(R, (1 - R)*active_parts)
z <- sample(0:(length(pi) - 1), N, replace = TRUE, prob = pi)
X <- matrix(nrow = N, ncol = 2)
X[z == 0,] <- runif(sum(z == 0) * 2)

mu <- matrix(nrow = K, ncol = 2)
Sigma <- array(dim = c(K,2,2))
L <- array(dim = c(K,2,2))

for(k in 1:K){
  mu[k,] <- runif(2)
  Sigma[k,,] <- (1/rgamma(1, 3, 2e-4))*diag(2)
  L[k,,] <- t(chol(Sigma[k,,]))
  X[z == k,] <- t(mu[k,] + L[k,,]%*%matrix(rnorm(sum(z == k)*2), nrow = 2))
}

sel_in_bounds <- X[,1] > 0 & X[,1] < 1 & X[,2] > 0 & X[,2] < 1
N <- sum(sel_in_bounds)
X <- X[sel_in_bounds,]
z <- z[sel_in_bounds]

data <- data.frame(x = X[,1], y = X[,2], z = factor(z))

png("output/sky_survey_analysis/simstudy_data.png", 
    width = 5, height = 5, units = 'in', res = 300)
ggplot() +
  geom_point(data = data %>% filter(z == '0'),
             aes(x = x, y = y), color = 'grey', alpha = 0.5, size = 0.1) +
  geom_point(data = data %>% filter(z != '0'),
             aes(x = x, y = y, color = z), size = 0.1) +
  geom_point(aes(x = mu[,1], y = mu[,2]),
             size = 2, color = 'black', shape = 4) +
  guides(color = 'none') +
  labs(title = "Simulated Sky Survey Data",
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
pb.close()
plot_grid$f_pe <- plot_grid$f_pe / nsamps

plot_grid$f_true <- R
for(k in 1:K){
  plot_grid$f_true <- plot_grid$f_true + (1 - R)*active_parts[k]*
    dnorm(plot_grid$x, mean = mu[k,1], sd = L[k,1,1])*
    dnorm(plot_grid$y, mean = mu[k,2], sd = L[k,2,2])
}

png("output/sky_survey_analysis/simstudy_randhist_density.png", 
    width = 10, height = 5, units = 'in', res = 300)
shared_breaks <- quantile(log(c(plot_grid$f_pe, plot_grid$f_true)), 
                          prob = c(0,seq(0.875,1,0.025)), names = FALSE)
p1 <- ggplot() + 
  geom_contour_filled(aes(x = x, y = y, z = log(f_pe)), breaks = shared_breaks, data = plot_grid) + 
  geom_point(aes(x = mu[,1], y = mu[,2]), shape = 4, size = 2) + 
  guides(fill = 'none') + labs(title = "Log of Posterior Expected Density")
p2 <- ggplot() + 
  geom_contour_filled(aes(x = x, y = y, z = log(f_true)), breaks = shared_breaks, data = plot_grid) + 
  geom_point(aes(x = mu[,1], y = mu[,2]), shape = 4, size = 2) + 
  guides(fill = 'none') + labs(title = "Log of True Density")
grid.arrange(p1, p2, nrow = 1)
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
# Find the BAND clusters
# ----------------------------------------
# target quantile is scientifically motivated
target_quantile <- 0.9

# Heuristically determine an appropriate minPts parameter
approximate_number_of_clusters <- 100
points_per_cluster <- (1 - target_quantile)*nobs / 
  approximate_number_of_clusters
minPts <- round(points_per_cluster / 4)
print(sprintf("Hueristic Suggested minPts: %d", minPts))


# determine an appropriate split_error_prob parameter
sperp_options <- seq(0.01, 0.25, 0.02)
delta_options <- c()
for(i in 1:length(sperp_options))
  delta_options[i] <- compute_lambda_delta(X, 
                                  f_samps, 
                                  target_quantile, 
                                  minPts, 
                                  split_err_prob = sperp_options[i])[2]
png("output/sky_survey_analysis/simstudy_delta_vs_sperp.png", 
    width = 10, height = 10, units = 'in', res = 300)
plot(sperp_options, delta_options, 
     main = "Delta vs. Split Err. Prob.", 
     xlab = NA, ylab = NA, type = 'l')
dev.off()

# There is no elbow, so stick with the default
split_err_prob <- 0.05


# Try finding an appropriate minPts algorithmically
minPts_options <- seq(10, 100, 10)
band_metrics <- data.frame()
for(minPts in minPts_options){
  clustering_samps <- density_based_clusterer(X, 
                                              f_samps, 
                                              cut_quantile = target_quantile,
                                              minPts = minPts,
                                              split_err_prob = split_err_prob)
  
  # Compute a quick, approximate point estimate based on this minPts value
  sel_active_pe <- colMeans(clustering_samps != 0) >= 0.5
  tmp_point_estimate <- rep(0, nobs)
  tmp_point_estimate[sel_active_pe] <- X[sel_active_pe,] %>%
    dist() %>%
    hclust(method = "single") %>%
    cutree(h = attr(clustering_samps, 'delta'))
  
  # Evaluate the quality of the point estimate
  band_metrics <- rbind(band_metrics,
                          list(minPts = minPts,
                               sensitivity = sensitivity(X, tmp_point_estimate, mu),
                               specificity = specificity(X, tmp_point_estimate, mu)
                          ))  
}
# Set minPts based on the observed metrics
minPts <- minPts_options[which.max(band_metrics$sensitivity + 
                                     band_metrics$specificity)]

png("output/sky_survey_analysis/band_metrics.png", 
    width = 5, height = 5, units = 'in', res = 300)
band_metrics %>% 
  pivot_longer(c(sensitivity, specificity), 
               names_to = 'metric',
               values_to = 'score') %>%
  ggplot() + 
  geom_line(aes(x = minPts, y = score, color = metric)) + 
  geom_vline(xintercept = minPts, color = 'red') +
  ylim(0,1) + 
  labs(title = "BAND Performance")
dev.off()

# There's seemingly no improvement in BAND from scanning over minPts,
# so stick with the default from the heuristic
minPts <- round(points_per_cluster / 4)


clustering_samps <- density_based_clusterer(X, 
                                            f_samps, 
                                            cut_quantile = target_quantile,
                                            minPts = minPts,
                                            split_err_prob = split_err_prob)

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

png("output/sky_survey_analysis/band_best_performance.png", 
    width = 15, height = 5, units = 'in', res = 300)
enriched_data_lb <- enrich_small_clusters(data, 'lb', size_lb = 25)
p1 <- ggplot() + 
  geom_point(aes(x = x, y = y, color = factor(lb)), 
             size = 0.1, data = data %>% filter(lb != 0)) + 
  geom_point(aes(x = x, y = y), 
             size = 0.1, alpha = 0.1, color = 'grey', 
             data = data %>% filter(lb == 0)) +
  stat_ellipse(aes(x = x, y = y, group = factor(lb)), 
               data = enriched_data_lb %>% filter(lb != 0)) +
  geom_point(aes(x = X1, y = X2), shape = 4, size = 2, data = data.frame(mu)) + 
  guides(color = 'none') + 
  labs(title = "2.5%-ile Lower Bound",
       subtitle = sprintf("minPts: %d, spep: %.2f, sensitivity: %.2f, specificity: %.2f", minPts, split_err_prob, sensitivity(X, data$lb, mu), specificity(X, data$lb, mu)),
       x = NULL, y = NULL)

enriched_data_pe <- enrich_small_clusters(data, 'pe', size_lb = 25)
p2 <- ggplot() + 
  geom_point(aes(x = x, y = y, color = factor(pe)), size = 0.1, data = data %>% filter(pe != 0)) + 
  geom_point(aes(x = x, y = y), size = 0.1, alpha = 0.1, color = 'grey', data = data %>% filter(pe == 0)) +
  stat_ellipse(aes(x = x, y = y, group = factor(pe)), data = enriched_data_pe %>% filter(pe != 0)) +
  geom_point(aes(x = X1, y = X2), shape = 4, size = 2, data = data.frame(mu)) + 
  guides(color = 'none') + 
  labs(title = "Point Estimate",
       subtitle = sprintf("minPts: %d, spep: %.2f, sensitivity: %.2f, specificity: %.2f", minPts, split_err_prob, sensitivity(X, data$pe, mu), specificity(X, data$pe, mu)),
       x = NULL, y = NULL)

enriched_data_ub <- enrich_small_clusters(data, 'ub', size_lb = 25)
p3 <- ggplot() + 
  geom_point(aes(x = x, y = y, color = factor(ub)), size = 0.1, data = data %>% filter(ub != 0)) + 
  geom_point(aes(x = x, y = y), size = 0.1, alpha = 0.1, color = 'grey', data = data %>% filter(ub == 0)) +
  stat_ellipse(aes(x = x, y = y, group = factor(ub)), data = enriched_data_ub %>% filter(ub != 0)) +
  geom_point(aes(x = X1, y = X2), shape = 4, size = 2, data = data.frame(mu)) + 
  guides(color = 'none') + 
  labs(title = "97.5%-ile Upper Bound",
       subtitle = sprintf("minPts: %d, spep: %.2f, sensitivity: %.2f, specificity: %.2f", minPts, split_err_prob, sensitivity(X, data$ub, mu), specificity(X, data$ub, mu)),
       x = NULL, y = NULL)
grid.arrange(p1, p2, p3, nrow = 1)
dev.off()


# ----------------------------------------
# Find the DBSCAN clusters
# ----------------------------------------
minPts_band <- minPts # save the value of minPts we were using for BAND
napprox <- 1000 # find epsilon based on an approximation of the full data
minPts_options <- seq(10, 100, 10)
tX <- t(X) # convenient to have this for operation broadcasting computations 

dbscan_metrics <- data.frame()
for(minPts in minPts_options){
  # Given minPts, determine an appropriate epsilon based on the target level
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
  tmp_point_estimate <- dbfit$cluster
  
  dbscan_metrics <- rbind(dbscan_metrics,
                          list(minPts = minPts,
                               sensitivity = sensitivity(X, tmp_point_estimate, mu),
                               specificity = specificity(X, tmp_point_estimate, mu)
                          ))
}
# Pick the best minPts based on dbscan_metrics
minPts <- minPts_options[which.max(dbscan_metrics$sensitivity + 
                                     dbscan_metrics$specificity)]
png("output/sky_survey_analysis/dbscan_metrics.png", 
    width = 5, height = 5, units = 'in', res = 300)
dbscan_metrics %>% 
  pivot_longer(c(sensitivity, specificity), 
               names_to = 'metric',
               values_to = 'score') %>%
  ggplot() + 
  geom_line(aes(x = minPts, y = score, color = metric)) + 
  geom_vline(xintercept = minPts, color = 'red') +
  ylim(0,1) + 
  labs(title = "DBSCAN Performance")
dev.off()


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
data$dbscan_best <- dbfit$cluster

dbfit <- dbscan(X, eps = eps, minPts = minPts_band, borderPoints = FALSE)
data$dbscan_bmp <- dbfit$cluster

png("output/sky_survey_analysis/dbscan_best_performance.png", 
    width = 5, height = 5, units = 'in', res = 300)
enriched_data <- enrich_small_clusters(data, 'dbscan_best', size_lb = 25)
ggplot() + 
  geom_point(aes(x = x, y = y, color = factor(dbscan_best)), 
             size = 0.1, data = data %>% filter(dbscan_best != 0)) + 
  geom_point(aes(x = x, y = y), 
             size = 0.1, alpha = 0.1, color = 'grey', 
             data = data %>% filter(dbscan_best == 0)) +
  stat_ellipse(aes(x = x, y = y, group = factor(dbscan_best)), 
               data = enriched_data %>% filter(dbscan_best != 0)) +
  geom_point(aes(x = X1, y = X2), shape = 4, size = 2, data = data.frame(mu)) + 
  guides(color = 'none') + 
  labs(title = "DBSCAN Estimated Clusters and True Locations",
       subtitle = sprintf("minPts: %d, eps: %.2f, sensitivity: %.2f, specificity: %.2f", minPts, eps, sensitivity(X, data$dbscan, mu), specificity(X, data$dbscan, mu)),
       x = NULL, y = NULL)
dev.off()

png("output/sky_survey_analysis/dbscan_band_minpts.png", 
    width = 5, height = 5, units = 'in', res = 300)
enriched_data <- enrich_small_clusters(data, 'dbscan_bmp', size_lb = 25)
ggplot() + 
  geom_point(aes(x = x, y = y, color = factor(dbscan_bmp)), 
             size = 0.1, data = data %>% filter(dbscan_bmp != 0)) + 
  geom_point(aes(x = x, y = y), 
             size = 0.1, alpha = 0.1, color = 'grey', 
             data = data %>% filter(dbscan_bmp == 0)) +
  stat_ellipse(aes(x = x, y = y, group = factor(dbscan_bmp)), 
               data = enriched_data %>% filter(dbscan_bmp != 0)) +
  geom_point(aes(x = X1, y = X2), shape = 4, size = 2, data = data.frame(mu)) + 
  guides(color = 'none') + 
  labs(title = "DBSCAN Estimated Clusters and True Locations",
       subtitle = sprintf("minPts: %d, eps: %.2f, sensitivity: %.2f, specificity: %.2f", minPts_band, eps, sensitivity(X, data$dbscan, mu), specificity(X, data$dbscan, mu)),
       x = NULL, y = NULL)
dev.off()

