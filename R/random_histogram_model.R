rbreaks <- function(range_x, res, alpha = 5){
  probs <- rgamma(res, alpha)
  probs <- c(0, cumsum(probs / sum(probs)))
  break_vals <- range_x[1] + diff(range_x)*probs
  break_vals[1] <- break_vals[1] * (1 - 1e-15)
  break_vals[res + 1] <- break_vals[res + 1] * (1 + 1e-15)
  break_vals
}
comp_area <- function(box){
  (box['xu'] - box['xl'])*(box['yu'] - box['yl'])
}
get_hist_data <- function(data, res, nhists){
  range_x <- with(data, range(x))
  range_y <- with(data, range(y))
  total_area <- diff(range_x)*diff(range_y)
  
  nobs <- nrow(data)
  # obs_alloc <- sample(1:nhists, nobs, replace = TRUE)
  breaks_x <- matrix(nrow = res + 1, ncol = nhists)
  breaks_y <- matrix(nrow = res + 1, ncol = nhists)
  box_areas <- matrix(nrow = res^2, ncol = nhists)
  box_counts <- matrix(nrow = res^2, ncol = nhists)
  for(h in 1:nhists){
    breaks_x[,h] <- with(data, rbreaks(range_x, res))
    breaks_y[,h] <- with(data, rbreaks(range_y, res))
    boxes <- cbind(expand.grid(xl = breaks_x[1:res,h],
                               yl = breaks_y[1:res,h]),
                   expand.grid(xu = breaks_x[2:(res + 1),h],
                               yu = breaks_y[2:(res + 1),h]))
    
    box_areas[,h] <- apply(boxes, 1, comp_area)
    # box_counts[,h] <- as.numeric(with(data[obs_alloc == h,], 
    #                                   table(cut(x, breaks_x[,h]), 
    #                                         cut(y, breaks_y[,h]))))
    box_counts[,h] <- as.numeric(with(data, 
                                      table(cut(x, breaks_x[,h]), 
                                            cut(y, breaks_y[,h]))))
  }
  return(list(breaks_x = breaks_x, 
              breaks_y = breaks_y, 
              box_areas = box_areas, 
              box_counts = box_counts))
}
sample_density <- function(targetPts,prior_alpha,hist_data){
  with(hist_data,
       {
         total_area <- diff(range(breaks_x))*
           diff(range(breaks_y))
         nboxes <- nrow(box_counts)
         nhists <- ncol(box_counts)
         res <- sqrt(nboxes)
         
         f_samp <- rep(0, nrow(targetPts))
         for(h in 1:nhists){
           box_probs <- rgamma(nboxes, box_counts[,h]/nhists + 
                                 prior_alpha*(box_areas[,h] / total_area))
           box_probs <- box_probs / sum(box_probs)
           box_density <- box_probs / box_areas[,h]
           
           x_ind <- with(targetPts, cut(x, breaks_x[,h], labels = FALSE))
           y_ind <- with(targetPts, cut(y, breaks_y[,h], labels = FALSE))
           
           kron_ind <- (y_ind - 1)*res + x_ind
           f_samp <- f_samp + box_density[kron_ind]
         }
         
         f_samp / nhists
       })
}