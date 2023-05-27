simplify_labels <- function(labels, threshold = 5){
  distinct <- unique(labels)
  for(didx in 1:length(distinct)){
    dl <- distinct[didx]
    dl_select <- labels == dl
    dl_count <- sum(dl_select)
    if(dl_count < threshold)
      labels[dl_select] <- 0
  }
  return(labels)
}

rearrange_labels <- function(labels){
  new_labels <- rep(0, length(labels))
  
  non0_labs <- which(labels != 0)
  labels <- labels[non0_labs]
  label_ranks <- rank(1/table(labels), ties.method = 'first')
  key <- as.numeric(names(label_ranks))
  value <- c(label_ranks)
  new_labels[non0_labs] <- sapply(labels,\(l) value[which(key == l)])
  return(new_labels)
}