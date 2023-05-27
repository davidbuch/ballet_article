splitrange <- function(v, n = 10){
  rv <- range(v)
  rv[1] + ((rv[2] - rv[1])/n)*c(-1e-16,1:n)
}

pt_prob <- function(bin_ids, 
                    which_bins = 1:nrow(bin_ids), 
                    depth = 1, 
                    max.depth = 7){
  bin_vec <- bin_ids[which_bins, depth]
  bin_counts <- table(bin_vec)
  
  # create a table of sampled probabilities
  probs <- rgamma(length(bin_counts), bin_counts + depth^2)
  probs <- probs / sum(probs)
  names(probs) <- names(bin_counts)
  
  prob_vals <- probs[as.character(bin_vec)]
  if(depth != max.depth){
    bins_to_split <- unique(bin_vec)
    for(bin in bins_to_split){
      bin_selector <- bin_vec == bin
      prob_vals[bin_selector] <- prob_vals[bin_selector]*
        Recall(bin_ids,which_bins[bin_selector], depth = depth + 1)
    }
  }
  return(prob_vals)
}

fit_tree_dirichlet <- function(X, S = 1000, tree_depth = 7){
  data <- data.frame(x = X[,1], y = X[,2])
  N <- nrow(data)
  
  bin_ids <- matrix(nrow = N, ncol = tree_depth)
  for(depth in 1:tree_depth){
    data <- data %>% 
      mutate(xbin = cut(x, splitrange(x, n = 2^depth)),
             ybin = cut(y, splitrange(y, n = 2^depth)),
             "xybin.{depth}" := factor(paste(xbin,ybin, sep = ", ")))
    bin_ids[,depth] <- as.integer(data[[paste0('xybin.',depth)]])
  }
  
  probs <- t(replicate(S, pt_prob(bin_ids)))
  return(probs)
}



