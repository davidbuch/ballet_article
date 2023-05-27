fn_from_apt <- function(res, xvals, max.resol){
  nobs <- nrow(xvals)
  obsdim <- ncol(xvals)
  nsims <- length(res$part_points_post_samples)
  
  fn_samps <- matrix(nrow = nsims, ncol = nobs)
  for(s in 1:nsims){
    if(s %% 10 == 0){
      print(s)
    }
    part.post.sample <- res$part_points_post_samples[[s]]
    terminal.part <- part.post.sample[which(part.post.sample[,"nu"] == Inf),]
    box_dens <- exp(terminal.part[,"logp"])*2^(terminal.part[,"level"])
    box_bounds <- terminal.part[,1:(2*obsdim),drop=FALSE]
    box_bounds[,2*(1:obsdim)] <- box_bounds[,2*(1:obsdim),drop=FALSE] + 1
    box_bounds <- box_bounds / 2^max.resol 
    
    for(i in 1:nobs){
      # print(apply(box_bounds,2,range))
      # graphic.params=gpar(fill="red")
      # plot.new()
      # grid.rect((box_bounds[,1]+box_bounds[,2])/2,(box_bounds[,3]+box_bounds[,4])/2,
      #           width=box_bounds[,2] - box_bounds[,1],height=box_bounds[,4] - box_bounds[,3])
      # grid.rect(xvals[i,1], xvals[i,2], width = 0.01, height = 0.01, gp = graphic.params)
      box_loc <- apply(sweep(box_bounds[,2*(1:obsdim) - 1,drop=FALSE],2,xvals[i,,drop=FALSE],FUN = "<") & 
        sweep(box_bounds[,2*(1:obsdim),drop=FALSE],2,xvals[i,,drop=FALSE],FUN = ">="), 1, all)
      
      # grid.rect((box_bounds[box_loc,1]+box_bounds[box_loc,2])/2,(box_bounds[box_loc,3]+box_bounds[box_loc,4])/2,
      #           width=box_bounds[box_loc,2] - box_bounds[box_loc,1],height=box_bounds[box_loc,4] - box_bounds[box_loc,3], gp = graphic.params)
      stopifnot(sum(box_loc) == 1)
      fn_samps[s,i] <- box_dens[box_loc]
    }
  }
  return(fn_samps)
}