library(RANN)
library(matrixStats)

density_based_clusterer <- function(x, 
                                    fsamps,
                                    cut_quantile = 0.2, 
                                    minPts = ceiling(log2(nrow(x))),
                                    split_err_prob = 0.01,
                                    lambda_delta,
                                    plot_dendro = FALSE)
{
  N <- ifelse(is.null(dim(x)), length(x), dim(x)[1])
  D <- ifelse(is.null(dim(x)), 1, dim(x)[2])
  stopifnot(ncol(fsamps) == N)
  S <- nrow(fsamps)
  
  if(missing(lambda_delta)){
    lambda_delta <- compute_lambda_delta(x,
                                         fsamps,
                                         cut_quantile = cut_quantile, 
                                         minPts = minPts,
                                         split_err_prob = split_err_prob)
    cat(sprintf("Lambda: %.2f, Delta: %.2f\n", 
                lambda_delta[1], lambda_delta[2]))
  }
  lambda <- lambda_delta[1]
  delta <- lambda_delta[2]
  
  labels <- matrix(nrow = S, ncol = N)
  pb <- txtProgressBar(width = 30, style = 3)
  for(s in 1:S){
    sel_core <- fsamps[s,] > lambda
    labels[s,] <- rep(0, N)
    # benchmark study suggests it's faster to recompute the distance matrix
    # than to save it and worry about as.matrix, as.dist function calls
    labels[s, sel_core] <- x[sel_core,] %>%
      dist() %>%
      hclust(method = "single") %>%
      cutree(h = delta)
    
    if(plot_dendro && (s %% round(S/5) == 0)){
      x[sel_core,] %>%
        dist() %>%
        hclust(method = "single") %>%
        plot(main = s)
      abline(h = delta, col = "red")
    }
    setTxtProgressBar(pb,value = s/S)
  }
  close(pb)
  
  attr(labels, 'lambda') <- lambda
  attr(labels, 'delta') <- delta
  
  labels
}

compute_lambda_delta <- function(x,
                                 fsamps,
                                 cut_quantile = 0.2, 
                                 minPts = ceiling(log2(nrow(x))),
                                 split_err_prob = 0.01)
{
  N <- ifelse(is.null(dim(x)), length(x), dim(x)[1])
  D <- ifelse(is.null(dim(x)), 1, dim(x)[2])
  stopifnot(ncol(fsamps) == N)

  Ef <- colMeans(fsamps)
  # set lambda here
  lambda <- quantile(Ef, probs = cut_quantile, names = FALSE)
  #lambda <- quantile(fsamps, probs = cut_quantile, names = FALSE)
  
  # # In the below, we select a radius (delta) for our topological tube so that, 
  # # if all core points had constant density greater than lambda within their 
  # # respective delta balls, then we would *usually find another point within 
  # # that point's delta ball*. The probability that there would be any core 
  # # points with empty 'level-lambda' delta balls is controlled by the 
  # # argument 'split_err_prob'.
  # 
  # Nc <- sum(colMeans(fsamps > lambda)) # estimate of number of core points
  # 
  # # Lower bound probability that an *individual* level-lambda delta ball will 
  # # have >=1 other core point
  # pr_gt1 <- (1 - split_err_prob)^1/Nc
  # 
  # # Upper bound probability that an individual level-lambda delta ball will have
  # # no other core points
  # pr_0 <- 1 - pr_gt1
  # 
  # # Upper bound probability that a specific 'other' core point will be outside 
  # # an individual level-lambda delta ball
  # p_miss <- pr_0^(1/(Nc - 1))
  # 
  # # Lower bound on the necessary probability coverage for an individual 
  # # level-lambda delta ball
  # p_cover <- 1 - p_miss
  # 
  # # Compute the necessary D-dimensional volume of the delta ball
  # V_D_delta <- p_cover/lambda
  # 
  # # Compute the radius of a d-dimensional ball with volume V_D_delta
  # delta <- ( (V_D_delta*gamma(D/2 + 1))^(1/D) ) / sqrt(3.14159)
  
  #expectedCorePts <- x[Ef > lambda,]
  #coreDists <- as.matrix(dist(expectedCorePts))
  #nnDists <- apply(coreDists, 1, \(r) sort(r)[minPts + 1])
  #delta <- quantile(nnDists, probs = 1 - split_err_prob, names = FALSE)
  
  delta <- delta_given_quantile(x, Ef, cut_quantile, minPts)
  
  return(c(lambda, delta))
}

# Our heuristic to compute delta given 
# - distX: the distance matrix of the entire dataset,
# - Ef: the estimated density of each point,
# and optional arguments:
# - cut_quantile: the fraction of noise points,
# - minPts: the size of the radius to ensure each active point has 
#           at least minPts many neighbors in the entire dataset.

delta_given_quantile <- function(x, Ef, cut_quantile=cut_quantile,
                                 minPts = ceiling(log2(nrow(x)))) {

  lambda <- quantile(Ef, probs = cut_quantile)
  k <- minPts + 1

  # Find the k-nearest neighbor distance for all elements in the active set.
  # The following code is equivalent to (but much faster):
  #  distX <- as.matrix(dist(x))
  #  nnDists <- apply(distX[Ef > lambda, ], 1, \(r) sort(r, partial=k)[k])
  nnDists <- rowMaxs(nn2(x, x[Ef > lambda, ], k, eps=1e-3, searchtype="priority")$nn.dists)
  
  delta  <- quantile(nnDists, probs = .99, names = FALSE)
  return(delta)
}
