library(here)
Rcpp::sourceCpp(here("src/ne_loss.cpp"))

compute_pst <- function(clustering_samps, 
                        color_map = function(x) ifelse(x == 0, 0, 1),
                        color_set = unique(color_map(0:ncol(clustering_samps)))){
  # Compute Pairwise Similarity Tensor
  # clustering_samps is an S by n matrix, where each of the S rows contains a 
  # length n allocation vector encoding a partition of [1:n]. 
  # color_map is a function which maps cluster components to their associated colors
  # by default, we make color_map a nonzero indicator function, to represent a baseline
  # vs. active part relationship.
  n <- ncol(clustering_samps)
  
  color_samps <- color_map(clustering_samps)
  pst <- array(0, dim = c(length(color_set), n, n))
  for(cidx in 1:length(color_set)){
    color <- color_set[cidx]
    for(i in 1:(n - 1)){
      for(j in (i + 1):n){
        pst[cidx, i, j] <- mean((clustering_samps[,i] == clustering_samps[,j])*
                (color_samps[,i] == color)*(color_samps[,j] == color))
      }
    }
    pst[cidx,,] <- pst[cidx,,] + t(pst[cidx,,])
  }
  return(pst)
}

compute_pdt <- function(clustering_samps, 
                        color_map = function(x) ifelse(x == 0, 0, 1),
                        color_set = unique(color_map(0:ncol(clustering_samps)))){
  # Compute Pairwise Similarity Tensor
  # clustering_samps is an S by n matrix, where each of the S rows contains a 
  # length n allocation vector encoding a partition of [1:n]. 
  # color_map is a function which maps cluster components to their associated colors
  # by default, we make color_map a nonzero indicator function, to represent a baseline
  # vs. active part relationship.
  n <- ncol(clustering_samps)
  
  color_samps <- color_map(clustering_samps)
  pdt <- array(0, dim = c(length(color_set), n, n))
  for(cidx in 1:length(color_set)){
    color <- color_set[cidx]
    for(i in 1:(n - 1)){
      for(j in (i + 1):n){
        pdt[cidx, i, j] <- mean((clustering_samps[,i] != clustering_samps[,j])*
                                  (color_samps[,i] == color)*(color_samps[,j] == color))
      }
    }
    pdt[cidx,,] <- pdt[cidx,,] + t(pdt[cidx,,])
  }
  return(pdt)
}

# ne_loss <- function(proposed_partition, 
#                  target_partition,
#                  color_map = function(x) ifelse(x == 0, 0, 1),
#                  color_set = unique(
#                    color_map(0:ncol(clustering_samps))
#                  ),
#                  sep_loss = 1, 
#                  join_loss = 1, 
#                  miscolor_loss = 1/4){
#   stopifnot(length(proposed_partition) == length(target_partition))
#   
#   g <- proposed_partition # cluster allocation vector `g`
#   gt <- target_partition
#   q <- color_map(g) # coloring vector `q`
#   qt <- color_map(gt)
#   
#   n <- length(g) # number of items to partition
#   
#   loss <- 0
#   for(i in 1:(n - 1)){
#     for(j in (i + 1):n){
#       use_binder <- (q[i] == qt[i]) && (q[j] == qt[j])
#       loss <- loss + 
#         miscolor_loss*((q[i] != qt[i]) + (q[j] != qt[j])) + 
#         sep_loss*(g[i] != g[j])*(gt[i] == gt[j])*use_binder + 
#         join_loss*(g[i] == g[j])*(gt[i] != gt[j])*use_binder
#     }
#   }
#   return(loss)
# }

ne_loss <- function(proposed_partition, 
                        target_partition,
                        color_map = function(x) ifelse(x == 0, 0, 1),
                        color_set = unique(
                          color_map(0:ncol(clustering_samps))
                        ),
                        sep_loss = 1, 
                        join_loss = 1, 
                        miscolor_loss = 1/4){
  stopifnot(length(proposed_partition) == length(target_partition))
  
  qprop <- color_map(proposed_partition) # coloring vector `q`
  qtarg <- color_map(target_partition)

  color_matched_indices <- which((qprop == qtarg) & qprop & qtarg)
  ne_loss_cpp(ghat = proposed_partition,
              g = target_partition, 
              qhat = qprop, 
              q = qtarg, 
              colormatch_indices = color_matched_indices, 
              sep_loss = sep_loss,
              join_loss = join_loss,
              miscolor_loss = miscolor_loss)
}

expected_loss <- function(proposed_partition,
                          clustering_samps,
                          color_map = function(x) ifelse(x == 0, 0, 1),
                          color_set = unique(
                            color_map(0:ncol(clustering_samps))
                          ),
                          sep_loss = 1, 
                          join_loss = 1, 
                          miscolor_loss = 1/4, 
                          quiet = TRUE){
  g <- proposed_partition # cluster allocation vector `g`
  q <- color_map(g) # coloring vector `q`
  n <- length(g) # number of items to partition
  
  color_samps <- color_map(clustering_samps)
  colors_match <- sweep(color_samps, 2, q, FUN = `==`)
  
  colors_match * clustering_samps
  
  prob_colors_not_matched <- colMeans(!colors_match)
  
  loss <- 0
  start_time <- Sys.time()
  for(i in 1:(n - 1)){
    for(j in (i + 1):n){
      use_binder <- colors_match[,i]*colors_match[,j]
      loss <- loss + 
        miscolor_loss*(prob_colors_not_matched[i] + prob_colors_not_matched[j]) + 
        sep_loss*(g[i] != g[j])*
          mean((clustering_samps[,i] == clustering_samps[,j])*use_binder) + 
        join_loss*(g[i] == g[j])*
          mean((clustering_samps[,i] != clustering_samps[,j])*use_binder)
    }
    if((i == 1) & !quiet){
      stop_time <- Sys.time()
      cat(sprintf("Computing the expected loss will take %.2f secs\n", 
              n * difftime(stop_time, start_time, units = 'secs')))
    }
  }
  return(loss)
}


expected_loss_marginal_form <- function(proposed_partition,
                                        clustering_samps,
                                        posterior_similarity_tensor,
                                        posterior_difference_tensor,
                                        color_map = function(x) ifelse(x == 0, 0, 1),
                                        color_set = unique(
                                          color_map(0:ncol(clustering_samps))
                                        ),
                                        sep_loss = 1, 
                                        join_loss = 1, 
                                        miscolor_loss = 1/4,
                                        always0_indcs = integer()){
  
  
  stopifnot(all(proposed_partition[always0_indcs] == 0))
  
  N <- length(proposed_partition)
  non0_indcs <- setdiff(1:N, always0_indcs)
  
  proposed_colors <- color_map(proposed_partition)
  color_samps <- color_map(clustering_samps)
  
  risk <- el_miscolor(proposed_colors, color_samps, miscolor_loss)

  risk <- risk + el_sepjoin(proposed_partition[non0_indcs],
                            proposed_colors[non0_indcs],
                            clustering_samps[,non0_indcs],
                            posterior_similarity_tensor,
                            posterior_difference_tensor,
                            color_set,
                            sep_loss,
                            join_loss)

  return(risk)
}

# Contribution of 'miscoloring' to the expected loss
el_miscolor <- function(proposed_colors,
                        color_samps,
                        miscolor_loss){
  n <- length(proposed_colors)

  prob_inactive <- colMeans(color_samps == 0)
  
  norm_risk <- sum(prob_inactive[proposed_colors == 1]) + 
    sum(1 - prob_inactive[proposed_colors == 0])
  
  return((n - 1)*miscolor_loss*norm_risk)
}

# Contribution of 'joins' and 'separations' within colors to the expected loss
el_sepjoin <- function(proposed_partition,
                       proposed_colors,
                       clustering_samps,
                       posterior_similarity_tensor,
                       posterior_difference_tensor,
                       color_set,
                       sep_loss,
                       join_loss){

  risk <- 0
  for(cidx in 1:length(color_set)){
    color <- color_set[cidx]
    Q_c <- which(proposed_colors == color)
    if(length(Q_c) >= 2){ 
      for(ii in 1:(length(Q_c) - 1)){
        for(jj in (ii + 1):length(Q_c)){
          i <- Q_c[ii]; j <- Q_c[jj]
          risk <- risk + 
            sep_loss*(proposed_partition[i] != proposed_partition[j])*
            posterior_similarity_tensor[cidx,i,j] + 
            join_loss*(proposed_partition[i] == proposed_partition[j])*
            posterior_difference_tensor[cidx,i,j]
        }
      }
    }
  }
  return(risk)
}


minimize_expected_loss <- function(clustering_samps, 
                                   color_map = function(x) ifelse(x == 0, 0, 1),
                                   color_set = unique(
                                     color_map(0:ncol(clustering_samps))
                                   ),
                                   sep_loss = 1,
                                   join_loss = 1,
                                   miscolor_loss = 1/4){
  # Compute the posterior expected loss for S distinct proposed partitions
  # e.g. from a multi-chain stochastic search
  # and return the proposed partition which has the minimum expected risk.
  n <- ncol(clustering_samps)
  color_samps <- color_map(clustering_samps)
  
  pst <- compute_pst(clustering_samps, color_map, color_set)
  pdt <- compute_pst(clustering_samps, color_map, color_set)
  
  S <- nrow(clustering_samps)
  sample_expected_risks <- vector(length = S)
  start_time <- Sys.time()
  for(s in 1:S){
    sample_expected_risks[s] <- expected_loss_marginal_form(clustering_samps[s,], 
                                                clustering_samps, 
                                                pst, pdt,
                                                color_map = color_map,
                                                color_set = color_set,
                                                sep_loss = sep_loss,
                                                join_loss = join_loss,
                                                miscolor_loss = miscolor_loss)
    if(s == 1){
      stop_time <- Sys.time()
      cat(sprintf("Minimizing the risk will take %.2f min\n", 
                  S * difftime(stop_time, start_time, units="mins")))
    }
  }
  
  cat(sprintf("Min Risk:  %f\n", min(sample_expected_risks)))
  
  return(clustering_samps[which.min(sample_expected_risks),])
}


