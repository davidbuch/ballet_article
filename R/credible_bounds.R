Rcpp::sourceCpp("src/ne_loss.cpp")
Rcpp::sourceCpp("src/subpartiton_min_max.cpp")

credible_bounds <- function(pe, 
                          clustering_samps,
                          prob = 0.95,
                          color_map = function(x) ifelse(x == 0, 0, 1),
                          color_set = unique(
                            color_map(0:ncol(clustering_samps))
                          ),
                          sep_loss = 1,
                          join_loss = 1,
                          miscolor_loss = 1/4){
  #' shell_width is the percentage up and down the distance quantile 
  #' we move to construct our shell
  S <- nrow(clustering_samps)
  N <- ncol(clustering_samps)
  
  color_pe <- color_map(pe)
  color_samps <- color_map(clustering_samps)
  
  clustering_distances <- vector(length = S)
  pb <- txtProgressBar(width = 30, style = 3)
  for(s in 1:S){
    colormatch_indices <- which((color_pe == color_samps[s,]) & 
                                  color_pe & color_samps[s,])
    
    clustering_distances[s] <- ne_loss_cpp(
      pe, 
      clustering_samps[s,],
      color_pe,
      color_samps[s,],
      colormatch_indices,
      sep_loss,
      join_loss,
      miscolor_loss
    )
    
    setTxtProgressBar(pb,value = s/S)
  }
  close(pb)
  

  epsilon_bound <- quantile(clustering_distances, prob, names = FALSE)
  ball_indcs <- which(clustering_distances <= epsilon_bound)
  
  credible_bounds <- list(
    min = min_subpartition(clustering_samps[ball_indcs,]),
    max = max_subpartition(clustering_samps[ball_indcs,])
  )
  
  return(credible_bounds)
}

# credible_bounds_active_inactive <- function(
#     x,
#     clustering_samps,
#     coverage_prob = 0.95,
#     color_map = function(x) ifelse(x == 0, 0, 1),
#     color_set = unique(
#       color_map(0:ncol(clustering_samps))
#     ),
#     sep_loss = 1,
#     join_loss = 1,
#     miscolor_loss = 1/4){
#   
#   N <- ncol(clustering_samps)
#   active_probabilities <- colMeans(clustering_samps != 0)
#   
#   sel_active_ub <- active_probabilities > (1 - coverage_prob) / 2
#   active_ub <- rep(0, N)
#   active_ub[sel_active_ub] <- x[sel_active_ub,] %>%
#     dist() %>%
#     hclust(method = "single") %>%
#     cutree(h = attr(clustering_samps, 'delta'))
#   
#   
#   sel_active_lb <- active_probabilities > 1 - ((1 - coverage_prob) / 2)
#   active_lb <- rep(0,N)
#   active_lb[sel_active_lb] <- x[sel_active_lb,] %>%
#     dist() %>%
#     hclust(method = "single") %>%
#     cutree(h = attr(clustering_samps, 'delta'))
#   
#   credible_bounds <- list(
#     lower = active_lb,
#     upper = active_ub
#   )
#   
#   return(credible_bounds)
# }


credible_ball_bounds_active_inactive <- function(
    x, 
    pe, 
    clustering_samps,
    prob = 0.95,
    color_map = function(x) ifelse(x == 0, 0, 1),
    color_set = unique(
      color_map(0:ncol(clustering_samps))
    ),
    sep_loss = 1,
    join_loss = 1,
    miscolor_loss = 1/4){
  
  nobs <- ncol(clustering_samps)
  nsamps <- nrow(clustering_samps)
  delta <- attr(clustering_samps,'delta')
  
  # Approximately compute the radius of the credible ball
  napprox <- 1000 # find epsilon
  row_samps <- sample(1:nsamps, napprox)
  partition_loss <- c()
  for(s in 1:napprox){
    partition_loss[s] <- ne_loss(pe, clustering_samps[row_samps[s],])
  }
  eps_star <- quantile(partition_loss, prob)
  
  # Now we will move 'up' and 'down' the lattice (by activity) to the bounds
  active_probabilities <- colMeans(clustering_samps != 0)
  
  middle_part <- partition_at_level(0.5, x, active_probabilities, delta)
  if(ne_loss(middle_part, pe) > eps_star)
    stop("Failed to find appropriate bounds - perhaps there is too little uncertainty?")

  # write a binary search tree to do this, I guess?
  
  # upper bound
  converged <- FALSE
  search_window <- c(0.5,1)
  while(!converged){
    test <- mean(search_window)
    part <- partition_at_level(test, x, active_probabilities, delta)
    if(ne_loss(part, pe) < eps_star){
      search_window[1] <- test
    }else{
      search_window[2] <- test
    }
    if(diff(search_window) < 0.001) converged <- TRUE
  }
  print(search_window)
  upper_bound <- partition_at_level(search_window[2], x, active_probabilities, delta)
  
  # lower bound
  converged <- FALSE
  search_window <- c(0,0.5)
  while(!converged){
    test <- mean(search_window)
    part <- partition_at_level(test, x, active_probabilities, delta)
    if(ne_loss(part, pe) > eps_star){
      search_window[1] <- test
    }else{
      search_window[2] <- test
    }
    if(diff(search_window) < 0.001) converged <- TRUE
  }
  print(search_window)
  lower_bound <- partition_at_level(search_window[1], x, active_probabilities, delta)
  
  # Note that we selected the *outer* representatives as bounds so they 
  # are conservative
  credible_bounds <- list(
    lower = lower_bound,
    upper = upper_bound
  )
  
  return(credible_bounds)
  
}

partition_at_level <- function(level,
                                   x,
                                   active_probabilities,
                                   delta){
  N <- nrow(x)
  sel <- active_probabilities > 1 - level
  
  partition <- rep(0, N)
  
  partition[sel] <- x[sel,] %>%
    dist() %>%
    hclust(method = "single") %>%
    cutree(h = delta)
  
  return(partition)
}
