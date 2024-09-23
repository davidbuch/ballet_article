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

## Cluster using BALLET - DPMM

# Obtain samples the posterior density evaluated at all the observations.

x <- matrix(x, ncol = 1)
fsamps <- matrix(nrow = S, ncol = nrow(x))
for(s in 1:S) fsamps[s,] <- PosteriorFunction(dpfit, s)(x)

## How should we set the cut-off level in level-set clustering? 
##
## Working assumption: 
###   Our clusters are regions of "high density",
##    well-separated by regions of "low density".
##    Thus our clusters are described by connected components
##    of {f \geq \lambda} for some density level \lambda > 0.
##
## To make our density cut-off \lambda comparable across different methods 
## (and density estimators), we translate our cut-off level into an 
## equivalent fraction of noise points parameter denoted by \nu \in (0,1) in 
## the paper and `cut_quantile` in the code. 
## Another intuition here is that ranking observations based on their density
## seems to be an easier problem than accurately estimating the density at each 
## value.

## Okay: let us see how DBSCAN sets is parameter `Eps` given a fixed value of `MinPts`. 
## Indeed, DBSCAN corresponds to level-set clustering using a kernel density estimate 
## with the "box" kernel. A short simulation study of DBSCAN heuristic and its 
## sensitivity can be found in the paper accompanying the dbscan R package 
## [Link](https://doi.org/10.18637/jss.v091.i01).
## 
## Basically, the DBSCAN heuristic suggests choosing some reasonably small 
## value of MinPts (typically 2*dim or dim + 1 or 2*dim - 1). Next given MinPts,
## the parameter `Eps` is tuned by finding a "sudden drop" in the graph of kNN 
## distances at the data points and their corresponding rank values. The intuition 
## here is that, the "noise points" will have a high kNN distance values while the 
## non-noise points will have low values. Thus a "drop" can be seen when we make 
## the transition. Further, is often shown that that the selected 
## `fraction of noise points` does not vary much under this heuristic as k is 
## slightly varied.

## Let us see what fraction of noise points would DBSCAN produce in this scenario:
library(dbscan)
library(kneedle)

# Assuming MinPts = 3 find Eps using the heuristic.
dists <- kNNdist(x, k=3, all=FALSE)
rank_and_dist <- kneedle(rank(dists), dists, decreasing = FALSE, concave=TRUE, sensitivity=2)
qplot(rank(dists), dists) + 
  geom_vline(xintercept=rank_and_dist[1], color="red")
cut_quantile <- 1 - rank_and_dist[1]/length(dists) 
# Thus 1.95% of the observations are declared as noise..


# We can follow a similar strategy 

## First, we obtain a density estimate from our data.
Ef <- matrixStats::colMedians(fsamps)

## Now we wish to find a "sudden rise" in the density values.
## In this example, since we expect the value of the 
## threshold to be close to zero, we have transformed the counts with a 
## log to clearly see this value. 
##
## Now let's see what value we obtain using our estimated densities
rank_and_den <- kneedle(rank(Ef), Ef)
plot(rank(Ef), log(Ef), xlab="ranks", ylab="sorted log(density)")
abline(v=rank_and_den[1], col="red")
cut_quantile <- rank_and_den[1]/length(Ef)  #0.017
# Thus we declare 1.7% of the observations as noise..

## Another diagnostic is the number of clusters. 
## If we obtain a single cluster.
## clusters <- level_set_clusters(x, Ef, cut_quantiles=seq(0,0.9,length.out=10), prefix="q")
## clustree(clusters, prefix="q")

lambda_delta <- compute_lambda_delta(x, fsamps, 
                                     cut_quantile = cut_quantile)

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
