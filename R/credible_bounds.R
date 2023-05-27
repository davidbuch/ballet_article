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
