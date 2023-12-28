library(RANN)
library(matrixStats)

density_based_clusterer <- function(x, 
                                    fsamps,
                                    cut_quantile = 0.1, 
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

  Ef <- colMedians(fsamps)
  lambda <- quantile(Ef, probs = cut_quantile, names = FALSE)
  delta <- delta_given_quantile(x, Ef, cut_quantile, minPts)
  
  return(c(lambda, delta))
}

# Our heuristic to compute delta given 
# - x: the observed n x d data matrix,
# - Ef: the estimated density of each point,
# - cut_quantile: the fraction of noise points
# and optional arguments:
# - minPts: the size of the radius to ensure each active point has 
#           at least minPts neighbors in the entire dataset.
# - split_err_prob: the quantile to choose among the nearest 
#                   neighbor distance from the top. Defaults to the 
#                   top 1%

delta_given_quantile <- function(x, Ef, cut_quantile,
                                 minPts = ceiling(log2(nrow(x))),
                                 split_err_prob = 0.01) {

  lambda <- quantile(Ef, probs = cut_quantile)
  k <- minPts + 1
  active <- Ef >= lambda

  # Find the k-nearest neighbor distance for all elements in the active set.
  # The following code is equivalent to (but much faster than):
  #  dists <- as.matrix(dist(x))
  #  nnDists <- apply(dists[active, ], 1, \(r) sort(r, partial=k)[k])
  nnDists <- rowMaxs(nn2(x, x[active, , drop=FALSE], k, 
                      eps=1e-3, searchtype="priority")$nn.dists)
  
  delta  <- quantile(nnDists, probs = 1-split_err_prob, names = FALSE)
  return(delta)
}


# The same function as above,
# but takes a list of quantiles as input and 
# efficiently computes the level-sets.
delta_given_quantiles <- function(x, Ef, cut_quantiles,
                                 minPts = ceiling(log2(nrow(x))),
                                 split_err_prob = 0.01) {

  cut_quantiles <- sort(cut_quantiles)
  lambda.zero <- quantile(Ef, probs = cut_quantiles[1])
  k <- minPts + 1
  active.ub <- Ef >= lambda.zero

  # Find the k-nearest neighbor distance for all elements in the active set.
  # The following code is equivalent to (but much faster than):
  #  dists <- as.matrix(dist(x))
  #  nnDists <- apply(dists[active, ], 1, \(r) sort(r, partial=k)[k])
  nnDists <- rowMaxs(nn2(x, x[active.ub, , drop=FALSE], k, 
                      eps=1e-3, searchtype="priority")$nn.dists)
  
  deltas <- numeric(length(cut_quantiles))

  for(i in seq_along(cut_quantiles)) {
    lambda <- quantile(Ef, probs = cut_quantiles[i])
    active <- Ef[active.ub] >= lambda
    deltas[i]  <- quantile(nnDists[active], 
                       probs = 1-split_err_prob, names = FALSE)
  }

  return(deltas)
}
