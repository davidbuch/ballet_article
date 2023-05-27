sensitivity <- function(X, labels_as_factor, target_locs){
  K <- nrow(target_locs)
  non_noise_labels <- setdiff(unique(labels_as_factor),'0')
  
  mu_in_hull <- rep(0, K)
  for(k in 1:K){
    for(kpe in non_noise_labels){
      if(sum(labels_as_factor == kpe) < 3) next
      ch <- geometry::convhulln(X[labels_as_factor == kpe,,drop=FALSE])
      status <- geometry::inhulln(ch, target_locs[k,,drop=FALSE])
      if(status){
        mu_in_hull[k] <- 1
        # p1 <- ggplot() + 
        #   geom_point(aes(x=x,y=y,color=z_pe), 
        #              size = 0.1,
        #              data = plot_data %>% filter(z_pe == kpe)) + 
        #   geom_point(aes(x=mu[k,1], y=mu[k,2]),
        #              size = 2, shape = 4) + 
        #   lims(x=c(0,1),y=c(0,1))
        # print(p1)
        break
      }
    }
  }
  mean(mu_in_hull)
}

specificity <- function(X, labels_as_factor, target_locs){
  non_noise_labels <- setdiff(unique(labels_as_factor), '0')
  
  hull_has_mu <- rep(0, length(non_noise_labels))
  for(label_id in 1:length(non_noise_labels)){
    kpe <- non_noise_labels[label_id]
    if(sum(labels_as_factor == kpe) < 3) next
    ch <- geometry::convhulln(X[labels_as_factor == kpe,,drop=FALSE])
    status <- geometry::inhulln(ch, target_locs)
    hull_has_mu[label_id] <- any(status)
  }
  mean(hull_has_mu)
}
