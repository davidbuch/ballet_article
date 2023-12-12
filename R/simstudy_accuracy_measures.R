sensitivity <- function(X, labels, target_locs){
  K <- nrow(target_locs)
  non_noise_labels <- setdiff(unique(labels),0)
  
  target_in_some_hull <- rep(0, K)
  for(k in 1:K){
    for(kpe in non_noise_labels){
      if(sum(labels == kpe) < 3) next
      ch <- geometry::convhulln(X[labels == kpe,,drop=FALSE])
      status <- geometry::inhulln(ch, target_locs[k,,drop=FALSE])
      if(status){
        target_in_some_hull[k] <- 1
        break
      }
    }
  }
  mean(target_in_some_hull)
}

enrich_small_clusters <- function(df, label_col, size_lb = 20, ntimes = 50){
  stopifnot(is.numeric(df %>% pull(!!rlang::sym(label_col))))

  cluster_info <- df %>% 
    group_by(!!rlang::sym(label_col)) %>%
    dplyr::summarise(mx = mean(x),
              my = mean(y),
              sdx = sd(x),
              sdy = sd(y),
              n = n())
  
  small_clusters <- (cluster_info %>% 
    filter(n < size_lb))[[label_col]]
  
  within_group_sds <- cluster_info[
    (cluster_info[[label_col]] != 0) & 
      (cluster_info[['n']] > size_lb),] %>%
    drop_na() %>% 
    select(sdx,sdy) %>% 
    colMeans()
  
  new_df <- df %>% select(c(x,y,rlang::sym(label_col)))
  if(length(small_clusters) > 0){
    for(i in 1:length(small_clusters)){
      new_df <- new_df[
        new_df[[label_col]] != small_clusters[i],]
      
      center <- as.numeric(cluster_info[
        cluster_info[[label_col]] == small_clusters[i],
        c('mx','my')
      ])
      
      fake_x <- matrix(
        rnorm(n = ntimes*2, 
              mean = rep(center, each = ntimes), 
              sd = rep(within_group_sds, each = ntimes)), 
        ncol = 2)
      
      tmp <- data.frame(
        x = fake_x[,1],
        y = fake_x[,2])
      tmp[label_col] <- small_clusters[i]
      
      new_df <- rbind(new_df, tmp)
    }
  }
  return(new_df)
}

specificity <- function(X, labels, target_locs, size_lb = 20, ntimes = 50){
  df <- data.frame(x = X[,1], y = X[,2], labs = labels)
  cluster_info <- df %>% 
    group_by(labs) %>%
    dplyr::summarise(mx = mean(x),
              my = mean(y),
              sdx = sd(x),
              sdy = sd(y),
              n = n())
  small_clusters <- (cluster_info %>% 
                       filter(n < size_lb))[['labs']]
  within_group_sds <- cluster_info[
    (cluster_info[['labs']] != 0) & 
      (cluster_info[['n']] > size_lb),] %>%
    drop_na() %>% 
    select(sdx,sdy) %>% 
    colMeans()
  
  non_noise_labels <- setdiff(unique(labels), 0)
  hull_has_mu <- rep(0, length(non_noise_labels))
  for(label_id in 1:length(non_noise_labels)){
    kpe <- non_noise_labels[label_id]
    if(sum(labels == kpe) > size_lb){
      ch <- geometry::convhulln(X[labels == kpe,])
    }else{
      center <- colMeans(X[labels == kpe,,drop=FALSE])
      fake_x <- matrix(rnorm(n = ncol(X)*ntimes, 
                             mean = rep(center, each = ntimes), 
                             sd = rep(within_group_sds, each = ntimes)), 
                       ncol = 2)
      ch <- geometry::convhulln(fake_x)
    }
    status <- geometry::inhulln(ch, target_locs)
    hull_has_mu[label_id] <- sum(status) == 1
  }
  mean(hull_has_mu)
}
