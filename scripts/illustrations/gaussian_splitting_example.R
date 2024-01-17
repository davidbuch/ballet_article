library(dirichletprocess)
library(salso)
library(matrixStats)
source("R/density_clusterer.R")
source("R/ne_parts_pair_counting.R")
source("R/salso_custom.R")
source("R/fn_from_apt.R")
source("R/scale_box.R")
source("R/choose_lambda.R")

set.seed(1234)

dir.create('output/illustrations', recursive = TRUE, showWarnings = FALSE)


# MCMC iterations when fitting Dirichlet process models
S <- 1000 

# Mixture Parameters
locations <- c(-1, 5, 9)
K <- 2
weights <- rep(1/K, K)

# Helper Functions
mixture_density <- function(x){
  # Evaluate the true mixture density at location x
  sum(weights * c(dunif(x, locations[1], locations[2]), 
                  dnorm(x, locations[3])))
}
mixture_density <- Vectorize(mixture_density)

generate_samples <- function(n){
  # Generate n samples from the true mixture density
  ifelse(runif(n) < weights[1], 
         runif(n, locations[1], locations[2]),
         rnorm(n, locations[3]))
}

posterior_expectation_f <- function(x, dpfit){
  # Evaluate the posterior expected density function (under dpfit, an output
  # from the dirichletprocess package) at locations x
  rowMeans(sapply(1:S, \(s) PosteriorFunction(dpfit,s)(x)))
}

x_grid <- seq(-2,12,0.2)
fx <- mixture_density(x_grid)

## 100 Samples
N <- 100
x <- generate_samples(N)
dpmod <- DirichletProcessGaussian(x)
dpfit <- Fit(dpmod, S)

clustering_samps <- t(simplify2array(dpfit$labelsChain))
clustering_point_estimate <- salso::salso(clustering_samps)

png("output/illustrations/univariate_gauss_unif_n100.png", 
    width = 8, height = 8, units = 'in', res = 300)
par(cex.axis=2)
plot(x_grid, fx, type = "l",
     xlab = NA, ylab = NA,
     lwd = 3)
lines(x_grid, posterior_expectation_f(x_grid, dpfit),
      lty = "dashed", col = "red", lwd = 3)
points(x, posterior_expectation_f(x, dpfit), col = clustering_point_estimate,
       cex = 2, pch = 16)
dev.off()

## 1000 Samples
N <- 1000
x <- generate_samples(N)
dpmod <- DirichletProcessGaussian(x)
dpfit <- Fit(dpmod, S)

clustering_samps <- t(simplify2array(dpfit$labelsChain))
clustering_point_estimate <- salso::salso(clustering_samps)

png("output/illustrations/univariate_gauss_unif_n1000.png", 
    width = 8, height = 8, units = 'in', res = 300)
par(cex.axis = 2)
plot(x_grid, fx, type = "l",
     xlab = NA, ylab = NA,
     lwd = 3)
lines(x_grid, posterior_expectation_f(x_grid, dpfit), 
      lty = "dashed", col = "red", lwd = 3)
sel <- sample(N, 100)
points(x[sel], posterior_expectation_f(x[sel], dpfit), 
       col = c('firebrick2', 'darkgreen', 'dodgerblue3', 'cyan2')[clustering_point_estimate[sel]], cex = 2, pch = clustering_point_estimate[sel])
dev.off()

## Cluster using BaLlet - DPMM
x <- matrix(x, ncol = 1)
fsamps <- matrix(nrow = S, ncol = nrow(x))
for(s in 1:S) fsamps[s,] <- PosteriorFunction(dpfit, s)(x)

Ef <- matrixStats::colMedians(fsamps)
  
clusters <- level_set_clusters(x,Ef, cut_quantiles=seq(0,0.9,length.out=10), prefix="q")
clustree(clusters, prefix="q")

lambda_delta <- compute_lambda_delta(x, fsamps, 
                                     cut_quantile = 0.05)

clustering_samps <- density_based_clusterer(x, fsamps, 
                                      lambda_delta = lambda_delta)
clustering_point_estimate <- salso_custom(clustering_samps)

png("output/illustrations/univariate_gauss_unif_ballet_dpmm.png", 
    width = 8, height = 8, units = 'in', res = 300)
par(cex.axis = 2)
plot(x_grid, fx, type = "l",
     xlab = NA, ylab = NA,
     lwd = 3)
lines(x_grid, posterior_expectation_f(x_grid, dpfit), 
      lty = "dashed", col = "red", lwd = 3)
abline(h = lambda_delta[1])
points(x[sel,], posterior_expectation_f(x[sel,], dpfit), 
       col = c('grey50', 'firebrick2', 'darkgreen')[clustering_point_estimate[sel] + 1], pch = (clustering_point_estimate[sel]), cex = 2)
dev.off()
