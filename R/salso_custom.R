# adopt 'ordered pair' convention for colored labels, flatten when needed
mode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

get_potential_allocations <- function(labels, label_options){
  labels <- na.omit(labels)
  allocation_options <- data.frame(q = NULL, z = NULL)
  for(i in 1:nrow(label_options)){
    extant_labels_in_color <- labels$z[labels$q == label_options$q[i]]
    new_label_in_color <- 
      min(
        max(c(0,extant_labels_in_color), na.rm = TRUE) + 1,
        label_options$max_z[i]
      )
    z_options_for_q <- data.frame(q = label_options$q[i],
                                  z = unique(c(extant_labels_in_color, 
                                               new_label_in_color)))
    allocation_options <- rbind(allocation_options,
                                z_options_for_q)
  }
  return(allocation_options)
}

# specific to our active-inactive label case
flatten_labels <- function(labels){
  flat_labels <- c()
  for(i in 1:nrow(labels)){
    flat_labels[i] <- ifelse(labels$q[i] == 'A', labels$z[i], 0)
  }
  return(flat_labels)
}

greedy_sequential_cluster_allocation <- function(clustering_samps,
                                                 pst, pdt,
                                                 label_options = data.frame(q = factor(c('I','A')), max_z = c(1,Inf)),
                                                 sep_loss = 1,
                                                 join_loss = 1,
                                                 miscolor_loss = 1/4,
                                                 always0_indcs = integer(),
                                                 verbose = FALSE){
  # we assume clustering_samps is the full SxN matrix of samples
  # but (to accomodate large N), pst and pdt need not be full 2xNxN tensors
  # and we have to introduce this weird indexing convention to accomodate that
  
  N <- ncol(clustering_samps)
  non0_indcs <- setdiff(1:N, always0_indcs)
  N_0 <- length(always0_indcs)
  N_non0 <- length(non0_indcs)
  
  # randomize allocation order
  # the pairwise similarity objects will be indexed *relatively*
  perm_rel <- sample(1:length(non0_indcs), length(non0_indcs))
  perm_glo <- non0_indcs[perm_rel]
  
  labels <- data.frame(q = factor(rep(NA, N), 
                                  levels = label_options[['q']]), 
                       z = as.integer(rep(NA, N)))
  labels[always0_indcs,] <- list(q = 'I', z = 1)
  
  # To-Do: Generalize for Multicolor Applications 
  flat_labels <- flatten_labels(labels)
  prob_inactive <- colMeans(clustering_samps == 0)
  
  if(verbose) pb <- txtProgressBar(width = 30, style = 3)
  for(i in 1:N_non0){
    ig <- perm_glo[i]
    ir <- perm_rel[i]
    
    best_allocation <- NA
    
    if(i == 1){
      # assign to first active cluster with probability 1 - Pr(ghat[ig] == 0)
      if(prob_inactive[ig] > 0.5){
        labels[ig,] <- list(q = 'I', z = 1)
        flat_labels[ig] <- 0
      }else{
        labels[ig,] <- list(q = 'A', z = 1)
        flat_labels[ig] <- 1
      }
    }else{
      potential_allocations <- get_potential_allocations(labels,label_options)
      v1 <- pst[2,perm_rel[1:(i - 1)],perm_rel[i]]
      v2 <- pdt[2,perm_rel[1:(i - 1)],perm_rel[i]]

      min_observed_delta_risk <- Inf
      for(p in 1:nrow(potential_allocations)){
        # get the cluster proposal and flatten it
        proposal <- potential_allocations[p,]
        flat_labels[ig] <- ifelse(proposal$q == 'A',proposal$z, 0)

        # TO-DO: is this correct? It seems overly conservative in forcing obs into
        # the inactive cluster, though perhaps that is sensible, as we can't
        # miscluster observations within the noise cluster
        delta_risk <- NA
        if(flat_labels[ig] == 0){
          delta_risk <- miscolor_loss*(i + N_0 - 1)*(1 - prob_inactive[ig])
        }else{
          same_as_prop <- flat_labels[perm_glo[1:(i - 1)]] == flat_labels[ig]
          delta_risk <- miscolor_loss*(i + N_0 - 1)*prob_inactive[ig] +
            sep_loss*sum( (!same_as_prop) * v1) +
            join_loss*sum( same_as_prop * v2)
        }

        if(delta_risk < min_observed_delta_risk){
          best_allocation <- flat_labels[ig]
          min_observed_delta_risk <- delta_risk
        }
      }
      
      # update labels (both the 'flat' and 'structured' copy)
      flat_labels[ig] <- best_allocation
      if(best_allocation == 0){
        labels[ig,] <- list(q = 'I', z = 1)
      }else{
        labels[ig,] <- list(q = 'A', z = best_allocation)
      }
    }
    
    if(verbose) setTxtProgressBar(pb, value = i/N_non0)
  }
  if(verbose) close(pb)
  
  return(labels)
}

cluster_sweetening <- function(labels, 
                               clustering_samps,
                               pst, pdt,
                               label_options = data.frame(q = factor(c('I','A')), max_z = c(1,Inf)),
                               sep_loss = 1,
                               join_loss = 1,
                               miscolor_loss = 1/4,
                               always0_indcs = integer(),
                               verbose = FALSE){
  # we assume clustering_samps is the full SxN matrix of samples
  # but (to accomodate large N), pst and pdt need not be full 2xNxN tensors
  # and we have to introduce this weird indexing convention to accomodate that
  stopifnot(nrow(labels) == ncol(clustering_samps))
  N <- ncol(clustering_samps)
  non0_indcs <- setdiff(1:N, always0_indcs)
  N_0 <- length(always0_indcs)
  N_non0 <- length(non0_indcs)
  
  # randomize allocation order
  # the pairwise similarity objects will be indexed *relatively*
  perm_rel <- sample(1:length(non0_indcs), length(non0_indcs))
  perm_glo <- non0_indcs[perm_rel]
  
  # To-Do: Generalize for Multicolor Applications 
  flat_labels <- flatten_labels(labels)
  prob_inactive <- colMeans(clustering_samps == 0)
  
  if(verbose) pb <- txtProgressBar(width = 30, style = 3)
  for(i in 1:N_non0){
    ig <- perm_glo[i]
    otherg <- setdiff(perm_glo, ig)
    ir <- perm_rel[i]
    otherr <- setdiff(perm_rel, ir)
    
    best_allocation <- NA
    
    potential_allocations <- get_potential_allocations(labels,label_options)
    v1 <- pst[2, otherr, ir]
    v2 <- pdt[2, otherr, ir]
    
    min_observed_delta_risk <- Inf
    for(p in 1:nrow(potential_allocations)){
      # get the cluster proposal and flatten it
      proposal <- potential_allocations[p,]
      flat_labels[ig] <- ifelse(proposal$q == 'A',proposal$z, 0)
      
      # TO-DO: is this correct? It seems overly conservative in forcing obs into
      # the inactive cluster, though perhaps that is sensible, as we can't
      # miscluster observations within the noise cluster
      delta_risk <- NA
      if(flat_labels[ig] == 0){
        delta_risk <- miscolor_loss*(N - 2)*(1 - prob_inactive[ig])
      }else{
        same_as_prop <- flat_labels[otherg] == flat_labels[ig]
        delta_risk <- miscolor_loss*(N - 2)*prob_inactive[ig] +
          sep_loss*sum( (!same_as_prop) * v1) +
          join_loss*sum( same_as_prop * v2)
      }
      
      if(delta_risk < min_observed_delta_risk){
        best_allocation <- flat_labels[ig]
        min_observed_delta_risk <- delta_risk
      }
    }
    
    # update labels (both the 'flat' and 'structured' copy)
    flat_labels[ig] <- best_allocation
    if(best_allocation == 0){
      labels[ig,] <- list(q = 'I', z = 1)
    }else{
      labels[ig,] <- list(q = 'A', z = best_allocation)
    }
    
    if(verbose) setTxtProgressBar(pb, value = i/N_non0)
  }
  if(verbose) close(pb)
  
  return(labels)
}

zealous_updates <- function(labels, 
                            clustering_samps,
                            pst, pdt,
                            label_options = data.frame(q = factor(c('I','A')), max_z = c(1,Inf)),
                            sep_loss = 1,
                            join_loss = 1,
                            miscolor_loss = 1/4,
                            always0_indcs = integer(),
                            verbose = FALSE){
  # we assume clustering_samps is the full SxN matrix of samples
  # but (to accomodate large N), pst and pdt need not be full 2xNxN tensors
  # and we have to introduce this weird indexing convention to accomodate that
  stopifnot(nrow(labels) == ncol(clustering_samps))
  N <- ncol(clustering_samps)
  non0_indcs <- setdiff(1:N, always0_indcs)
  labels_update <- labels
  distinct_clusters <- distinct(labels)
  
  # get list of current clusters
  distinct_clusters <- distinct(labels)
  
  # randomly select a cluster to destroy
  target <- distinct_clusters %>% slice_sample(n = 1)

  # create flat versions of the labels and target
  flat_labels <- flatten_labels(labels)
  flat_labels_update <- flatten_labels(labels_update)
  flat_target <- flatten_labels(target)
  
  # get the indices which are assigned to that cluster to reallocate
  realloc <- setdiff(which(flat_labels == flat_target), always0_indcs)
  preserve <- setdiff(which(flat_labels != flat_target), always0_indcs)
  labels_update[realloc,] <- list(q = NA, z = NA)
  flat_labels_update[realloc] <- NA
  
  N_0 <- length(always0_indcs)
  N_p <- length(preserve)
  N_r <- length(realloc)
  
  # randomize allocation order
  # the pairwise similarity objects will be indexed *relatively*
  perm_rel <- sample(1:N_r, N_r)
  perm_glo <- realloc[perm_rel]
  preserve_rel <- match(preserve,non0_indcs)
  
  # pre-compute the probability of each point being assigned to 'inactive'
  prob_inactive <- colMeans(clustering_samps == 0)
  
  if(verbose) pb <- txtProgressBar(width = 30, style = 3)
  for(i in 1:N_r){
    ig <- perm_glo[i]
    otherg <- c(preserve, perm_glo[1:(i - 1)]) #setdiff(perm_glo, ig)
    ir <- perm_rel[i]
    otherr <- c(preserve_rel, perm_rel[1:(i - 1)]) # setdiff(perm_rel, ir)
    
    best_allocation <- NA
    
    potential_allocations <- get_potential_allocations(labels_update, label_options)
    v1 <- pst[2, otherr, ir]
    v2 <- pdt[2, otherr, ir]
    
    min_observed_delta_risk <- Inf
    for(p in 1:nrow(potential_allocations)){
      # get the cluster proposal and flatten it
      proposal <- potential_allocations[p,]
      flat_labels_update[ig] <- ifelse(proposal$q == 'A',proposal$z, 0)
      
      # TO-DO: is this correct? It seems overly conservative in forcing obs into
      # the inactive cluster, though perhaps that is sensible, as we can't
      # miscluster observations within the noise cluster
      delta_risk <- NA
      if(flat_labels_update[ig] == 0){
        delta_risk <- miscolor_loss*(i + N_0 + N_p - 1)*(1 - prob_inactive[ig])
      }else{
        same_as_prop <- flat_labels[otherg] == flat_labels[ig]
        delta_risk <- miscolor_loss*(i + N_0 + N_p - 1)*prob_inactive[ig] +
          sep_loss*sum( (!same_as_prop) * v1) +
          join_loss*sum( same_as_prop * v2)
      }
      
      if(delta_risk < min_observed_delta_risk){
        best_allocation <- flat_labels[ig]
        min_observed_delta_risk <- delta_risk
      }
    }
    
    # update labels (both the 'flat' and 'structured' copy)
    flat_labels_update[ig] <- best_allocation
    if(best_allocation == 0){
      labels_update[ig,] <- list(q = 'I', z = 1)
    }else{
      labels_update[ig,] <- list(q = 'A', z = best_allocation)
    }
    
    if(verbose) setTxtProgressBar(pb, value = i/N_r)
  }
  if(verbose) close(pb)
  
  if(verbose) cat("Comparing risk of previous clustering and zealous update...\n")
  risk_prev <- expected_loss_marginal_form(flatten_labels(labels),
                                     clustering_samps,
                                     pst,pdt,
                                     always0_indcs = always0_indcs)
  risk_updated <- expected_loss_marginal_form(flatten_labels(labels_update),
                                        clustering_samps,
                                        pst,pdt,
                                        always0_indcs = always0_indcs)
  if(risk_updated < risk_prev){
    if(verbose) cat("Updating labels!\n")
    labels <- labels_update
  }
  
  return(labels)
}

random_sequential_cluster_allocation <- function(n_labels,
                                                 label_options = data.frame(q = factor(c('I','A')), max_z = c(1,Inf)),
                                                 always0_indcs = integer(), 
                                                 verbose = FALSE){
  # we assume clustering_samps is the full SxN matrix of samples
  # but (to accomodate large N), pst and pdt need not be full 2xNxN tensors
  # and we have to introduce this weird indexing convention to accomodate that
  N <- n_labels
  non0_indcs <- setdiff(1:N, always0_indcs)
  N_0 <- length(always0_indcs)
  N_non0 <- length(non0_indcs)
  
  labels <- data.frame(q = factor(rep(NA, N), 
                                  levels = label_options[['q']]), 
                       z = as.integer(rep(NA, N)))
  labels[always0_indcs,] <- list(q = 'I', z = 1)
  
  # To-Do: Generalize for Multicolor Applications 
  flat_labels <- flatten_labels(labels)

  if(verbose) pb <- txtProgressBar(width = 30, style = 3)
  for(i in 1:N_non0){
    ig <- non0_indcs[i]
    
    potential_allocations <- get_potential_allocations(labels,label_options)
    p <- sample(1:nrow(potential_allocations), 1)
    proposal <- potential_allocations[p,]
    
    # update labels (both the 'flat' and 'structured' copy)
    random_allocation <- ifelse(proposal$q == 'A',proposal$z, 0)
    
    flat_labels[ig] <- random_allocation
    if(random_allocation == 0){
      labels[ig,] <- list(q = 'I', z = 1)
    }else{
      labels[ig,] <- list(q = 'A', z = random_allocation)
    }
    
    if(verbose) setTxtProgressBar(pb, value = i/N_non0)
  }
  if(verbose) close(pb)
  
  return(labels)
}

salso_single_run <- function(clustering_samps,
                             pst, pdt,
                             nZealous = 10,
                             label_options = data.frame(q = factor(c('I','A')), max_z = c(1,Inf)),
                             sep_loss = 1,
                             join_loss = 1,
                             miscolor_loss = 1/4,
                             always0_indcs = integer(),
                             verbose = FALSE){
  # initialize the cluster allocations
  if(runif(1) < 0.5){
    if(verbose) cat("Attempting random initialization... \n")
    labels <- random_sequential_cluster_allocation(n_labels = ncol(clustering_samps),
                                                   label_options = label_options,
                                                   always0_indcs = always0_indcs,
                                                   verbose = verbose)
  }else{
    if(verbose) cat("Attempting greedy initialization... \n")
    labels <- greedy_sequential_cluster_allocation(clustering_samps,
                                                   pst, pdt,
                                                   label_options = label_options,
                                                   sep_loss = sep_loss,
                                                   join_loss = join_loss,
                                                   miscolor_loss = miscolor_loss,
                                                   always0_indcs = always0_indcs,
                                                   verbose = verbose)
  }
  
  # 'sweeten' until no changes are observed
  sweetening_counter <- 0
  while(TRUE){
    if(verbose) cat("Sweetening... \n")
    labels_cache <- labels
    labels <- cluster_sweetening(labels, clustering_samps,
                                 pst, pdt,
                                 label_options = label_options,
                                 sep_loss = sep_loss,
                                 join_loss = join_loss,
                                 miscolor_loss = miscolor_loss,
                                 always0_indcs = always0_indcs,
                                 verbose = verbose)
    
    changed_labels <- sum(flatten_labels(labels) != flatten_labels(labels_cache))
    sweetening_counter <- sweetening_counter + 1
    if(verbose) 
      cat(sprintf("Completed %d round(s) of sweetening, changing %d labels! \n", sweetening_counter, changed_labels))
    if(changed_labels == 0) break
  }
  
  for(zcount in 1:nZealous){
    if(verbose) cat("Performing zealous update... \n")
    labels <- zealous_updates(labels, clustering_samps,
                              pst, pdt,
                              label_options = label_options,
                              sep_loss = sep_loss,
                              join_loss = join_loss,
                              miscolor_loss = miscolor_loss,
                              always0_indcs = always0_indcs,
                              verbose = verbose)
    if(verbose) cat(sprintf("Completed %d round(s) of zealous updates! \n", zcount))
  }
  
  return(labels)
}


salso_custom <- function(clustering_samps,
                  pst, pdt,
                  nRuns = 10,
                  nZealous = 10,
                  label_options = data.frame(q = factor(c('I','A')), max_z = c(1,Inf)),
                  sep_loss = 1,
                  join_loss = 1,
                  miscolor_loss = 1/4,
                  always0_indcs = integer(),
                  verbose = FALSE){
  if(missing(pst))
    pst <- compute_pst(clustering_samps)
  if(missing(pdt)) 
    pdt <- compute_pdt(clustering_samps)
  
  labels_list <- list()
  risk_list <- c()
  for(n in 1:nRuns){
    labels_list[[n]] <- salso_single_run(clustering_samps = clustering_samps,
                                         pst = pst, pdt = pdt,
                                         nZealous = nZealous,
                                         label_options = label_options,
                                         sep_loss = sep_loss,
                                         join_loss = join_loss,
                                         miscolor_loss = miscolor_loss,
                                         always0_indcs = always0_indcs,
                                         verbose = verbose)
    
    risk_list[n] <- expected_loss_marginal_form(
      flatten_labels(labels_list[[n]]),
      clustering_samps = clustering_samps,
      posterior_similarity_tensor = pst, posterior_difference_tensor = pdt,
      always0_indcs = always0_indcs
      )
  }
  if(verbose){
    cat("Risk Values: ")
    cat(risk_list)
    cat("\n")
  }
  return(flatten_labels(labels_list[[which.min(risk_list)]]))
}

