library(tidyverse)
source("R/random_histogram_model.R")
source("R/density_clusterer.R")
source("R/ne_parts_pair_counting.R")
source("R/salso_custom.R")
source("R/credible_bounds.R")
source("R/simstudy_accuracy_measures.R")

random_seed <- 1234
set.seed(random_seed)

N <- 40000 # n obs
R <- 0.90 # noise fraction
K <- 42 # number of non-noise clusters

# simulate active allocations from a dirichlet dsn
active_parts <- rgamma(K, 3)
active_parts <- active_parts / sum(active_parts)
barplot(active_parts)

#pi <- c(R, rep((1 - R) / K, K))
pi <- c(R, (1 - R)*active_parts)
z <- sample(0:(length(pi) - 1), N, replace = TRUE, prob = pi)
X <- matrix(nrow = N, ncol = 2)
X[z == 0,] <- runif(sum(z == 0) * 2)

mu <- matrix(nrow = K, ncol = 2)
Sigma <- array(dim = c(K,2,2))
L <- array(dim = c(K,2,2))

for(k in 1:K){
  mu[k,] <- runif(2)
  #Sigma[k,,] <- solve(rWishart(1,5,1e4*diag(2))[,,1])
  Sigma[k,,] <- (1/rgamma(1, 3, 2e-4))*diag(2)
  L[k,,] <- t(chol(Sigma[k,,]))
  X[z == k,] <- t(mu[k,] + L[k,,]%*%matrix(rnorm(sum(z == k)*2), nrow = 2))
}

sel_in_bounds <- X[,1] > 0 & X[,1] < 1 & X[,2] > 0 & X[,2] < 1
X <- X[sel_in_bounds,]
z <- z[sel_in_bounds]

data <- data.frame(x = X[,1], y = X[,2], z = factor(z))

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

## Fit Random Histogram
nobs <- nrow(data)
nsamps <- 1000

prior_alpha <- 1
res <- 50
nhists <- 50

hist_data <- get_hist_data(data, res, nhists)

# Sample from the posterior of the density function
f_samps <- matrix(nrow = nsamps, ncol = nrow(data))
pb <- txtProgressBar(width = 30, style = 3)
for(s in 1:nsamps){
  f_samps[s,] <- sample_density(data,
                                prior_alpha,
                                hist_data)
  setTxtProgressBar(pb, s/nsamps)
}
close(pb)

# Find the density based clusters
target_quantile <- 0.9
approximate_number_of_clusters <- 100
points_per_cluster <- (1 - target_quantile)*nobs / 
  approximate_number_of_clusters

minPts <- round(points_per_cluster / 2)
split_err_prob <- 0.05
clustering_samps <- density_based_clusterer(X, 
                                            f_samps, 
                                            cut_quantile = target_quantile,
                                            minPts = minPts,
                                            split_err_prob = split_err_prob)

always0_indcs <- which(apply(clustering_samps, 2, \(v) mean(v == 0) > 0.8))
non0_clustering_samps <- clustering_samps[,-always0_indcs]

posterior_similarity_tensor <- compute_pst(non0_clustering_samps)
posterior_difference_tensor <- compute_pdt(non0_clustering_samps)

data$pe <- factor(salso_custom(clustering_samps,
                        pst = posterior_similarity_tensor,
                        pdt = posterior_difference_tensor,
                        always0_indcs = always0_indcs))

ggplot() +
  geom_point(data = data %>% filter(pe == '0'),
             aes(x = x, y = y), color = 'grey', alpha = 0.5, size = 0.1) +
  geom_point(data = data %>% filter(pe != '0'),
             aes(x = x, y = y, color = pe), size = 0.1) +
  microViz::stat_chull(data = data %>% filter(pe != '0'),
                       aes(x = x, y = y, color = pe)) + 
  geom_point(aes(x = mu[,1], y = mu[,2]),
             size = 2, color = 'black', shape = 4) +
  guides(color = 'none') +
  labs(title = "BAND Estimated Clusters and True Locations",
       subtitle = sprintf("minPts = %d, split error prob. = %.2f", 
                          minPts, split_err_prob),
       x = NULL, y = NULL)

sensitivity(X, data$pe, mu)
specificity(X, data$pe, mu)

data$lb <- min_subpartition(clustering_samps)
sensitivity(X, data$lb, mu)
specificity(X, data$lb, mu)

# DBSCAN scored 0.85, 0.673



# Plot the Density Estimate on a Grid
grid_res <- 100
plot_grid <- with(data, 
                  expand.grid(x = seq(min(x), max(x)*(1-1e-8), length.out = grid_res),
                              y = seq(min(y), max(y)*(1-1e-8), length.out = grid_res))
)
plot_grid$f_pe <- rep(0, nrow(plot_grid))
for(s in 1:nsamps){
  plot_grid$f_pe <- plot_grid$f_pe + 
    sample_density(plot_grid,
                   prior_alpha,
                   hist_data
    )
}
plot_grid$f_pe <- plot_grid$f_pe / nsamps

ggplot() + 
  geom_contour_filled(aes(x = x, y = y, z = f_pe), data = plot_grid) + 
  geom_point(aes(x = mu[,1], y = mu[,2]), shape = 4, size = 2) + 
  guides(fill = 'none')













# find credible range (go back to unmasked labels)
# clustering_samps <- readRDS("fitted_models/simulated_data/galaxy_study/tree_dirichlet_clustering_samps.rds")
# labels_expanded <- readRDS("fitted_models/simulated_data/galaxy_study/tree_dirichlet_clustering_pe.rds")

lower_upper_bounds <- credible_range(labels_expanded,
                                     clustering_samps)

plot_data$z_lower <- factor(lower_upper_bounds[1,])
plot_data$z_upper <- factor(lower_upper_bounds[2,])

test_results <- matrix(nrow = 3, ncol = 2)
test_results[1,1] <- sensitivity(X, plot_data$z_pe, mu)
test_results[1,2] <- specificity(X, plot_data$z_pe, mu)

test_results[2,1] <- sensitivity(X, plot_data$z_lower, mu)
test_results[2,2] <- specificity(X, plot_data$z_lower, mu)

test_results[3,1] <- sensitivity(X, plot_data$z_upper, mu)
test_results[3,2] <- specificity(X, plot_data$z_upper, mu)

test_results <- data.frame(seed = random_seed,
                           type = c("PE",
                                    "Lower",
                                    "Upper"),
                           sensitivity = test_results[,1],
                           specificity = test_results[,2])

dir_name <- "fitted_models/simulated_data/galaxy_study"
if(!dir.exists(dir_name)) dir.create(dir_name)

filename <- sprintf(paste0(dir_name, "/accuracy_%d.rds"), random_seed)
saveRDS(test_results, filename)



## DBSCAN
library(dbscan)
row_samps <- sample(1:nobs, 1000)

dbscan_metrics <- data.frame()
for(minPts in seq(10, 100, 10)){
tX <- t(X)
minPts_NN_dists <- c()
for(ridx in row_samps){
  x <- tX[,ridx]
  dists <- sqrt(colSums((tX - x)^2))
  minPts_NN_dists <- c(minPts_NN_dists,
                       sort(dists, partial = minPts)[minPts])
}
minPts_NN_dists <- sort(minPts_NN_dists)
#plot(1:1000, minPts_NN_dists, type = "l")
eps <- minPts_NN_dists[round(1000*(1 - target_quantile))]
#abline(h = eps, col = 'red')

dbfit <- dbscan(X, eps = eps, minPts = minPts, borderPoints = FALSE)
data$dbscan <- dbfit$cluster

# ggplot() +
#   geom_point(data = data %>% filter(dbscan == '0'),
#              aes(x = x, y = y), color = 'grey', alpha = 0.5, size = 0.1) +
#   geom_point(data = data %>% filter(dbscan != '0'),
#              aes(x = x, y = y, color = factor(dbscan)), size = 0.1) + 
#   geom_point(aes(x=mu[,1], y=mu[,2]),
#              size = 2, shape = 4) +
#   guides(color = 'none') +
#   labs(title = "DBSCAN Estimated Clusters and True Locations",
#        subtitle = sprintf("minPts = %d, eps = %.2f", minPts, eps),
#        x = NULL, y = NULL)


dbscan_metrics <- rbind(dbscan_metrics,
                        list(minPts = minPts,
                             sensitivity = sensitivity(X, data$dbscan, mu),
                             specificity = specificity(X, data$dbscan, mu)
))
}
dbscan_metrics %>% 
  pivot_longer(c(sensitivity, specificity), 
               names_to = 'metric',
               values_to = 'score') %>%
  ggplot() + 
  geom_line(aes(x = minPts, y = score, color = metric)) + 
  labs(title = "DBSCAN Performance")
