library(tidyverse)
library(dbscan)
source("R/random_histogram_model.R")
source("R/density_clusterer.R")
source("R/ne_parts_pair_counting.R")
source("R/salso_custom.R")
source("R/credible_bounds.R")
source("R/simstudy_accuracy_measures.R")

random_seed <- as.integer(Sys.getenv(SLURM_ARRAY_TASK_ID))
set.seed(random_seed)

# ----------------------------------------
# Simulate a New Dataset
# ----------------------------------------
N <- 40000 # n obs
R <- 0.90 # noise fraction
K <- 42 # number of non-noise clusters

# simulate active allocations from a dirichlet dsn
active_parts <- rgamma(K, 3)
active_parts <- active_parts / sum(active_parts)

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

# ----------------------------------------
# Fit the Random Histogram Model
# ----------------------------------------
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
  tryCatch({
    f_samps[s,] <- sample_density(data,
                                  prior_alpha,
                                  hist_data)
  }, error = function(e){s <- s - 1})
  setTxtProgressBar(pb, s/nsamps)
}
close(pb)

# ---------------------------------------------------------
# Find the BAND clusters (Point Estimate and Bounds)
# ---------------------------------------------------------
# Find the density based clusters
target_quantile <- 0.9
minPts <- 20 # determined from pre-simulation
split_err_prob <- 0.05 # determined from pre-simulation

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

bounds <- credible_bounds_active_inactive(X, clustering_samps)
data$lb <- bounds$lower
data$ub <- bounds$upper

# ----------------------------------------
# Find the DBSCAN clusters
# ----------------------------------------
napprox <- 1000 # find epsilon based on an approximation of the full data
tX <- t(X) # convenient to have this for operation broadcasting computations 
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

# ----------------------------------------
# Save the Results
# ----------------------------------------
test_results <- data.frame(
  seed = random_seed,
  type = c("DBSCAN",
           "BAND_LB",
           "BAND_PE",
           "BAND_UB"
  ),
  sensitivity = c(
    sensitivity(X, data$dbscan, mu),
    sensitivity(X, data$lb, mu),
    sensitivity(X, data$pe, mu),
    sensitivity(X, data$ub, mu)
  ),
  specificity = c(
    specificity(X, data$dbscan, mu),
    specificity(X, data$lb, mu),
    specificity(X, data$pe, mu),
    specificity(X, data$ub, mu)
  )
)

dir_name <- "output/sky_survey_analysis/sim_study"
dir.create(dir_name, showWarnings = FALSE, recursive = TRUE)
filename <- sprintf(paste0(dir_name, "/accuracy_%d.rds"), random_seed)
saveRDS(test_results, filename)