# A collection of functions to plot an estimate of the 
# cluster tree across various levels of lambda.

source("R/density_clusterer.R")
library(usedist)
library(tidygraph)

 
# This function estimates the level set clustering
# based on a density estimate fest (length n vector),  
# a distance matrix distX (n x n matrix), and 
# cut_quantile (a number between 0 and 1).
level_set_clust <- function(cut_quantile, distX, fest, delta=NULL) {
    distX <- as.matrix(distX)
    
    if(is.null(delta)) {
      delta <- delta_given_quantile(distX, 
                                    fest, 
                                    cut_quantile=cut_quantile)
    }
  
    labels <- rep(0, nrow(distX))
    active <- fest >= quantile(fest, cut_quantile)
    
    if(sum(active) >= 2) {
      # If there are atleast two points in the active set,
      # cluster them using single linkage clustering, at 
      # threshold delta.
      distX |>
        usedist::dist_subset(active) |>
          hclust(method = "single") |> 
          cutree(h = delta) -> labels[active]
    } else if (sum(active) == 1) {
      labels[active] <- 1
    }
    return(labels)
} 

# This function estimates a cluster tree based on 
# data X (n x p matrix), a density estimate fest (length n vector),
# and a collection of quantiles cut_quantiles.
level_set_clusters <- function(X, fest, 
                            cut_quantiles=seq(0, 1, by = 0.1), 
                           delta=NULL,
                           prefix="q") {
    n <- nrow(X)
    distX <- as.matrix(dist(X))
    
    clusters <- vapply(X=cut_quantiles, FUN=level_set_clust, 
                       FUN.VALUE = rep(0, n), 
                       distX, fest, delta)
    colnames(clusters) <- sprintf("%s%.3f", prefix, cut_quantiles)
    clusters
}


#' Get tree nodes
#'
#' Extract the nodes from a set of clusterings and add relevant attributes
#'
#' @param clusterings numeric matrix containing clustering information, each
#' column contains clustering at a separate resolution
#' @param metadata data.frame containing metadata on each sample that can be
#' used as node aesthetics
#' @param prefix string indicating columns containing clustering information
#' @param node_aes_list nested list containing node aesthetics
#'
#' @keywords internal
#'
#' @return data.frame containing node information
get_tree_nodes_p <- function(clusterings, prefix, metadata, node_aes_list) {
  
  nodes <- lapply(colnames(clusterings), function(res) {
    clustering <- clusterings[, res]
    
    clusters <- sort(unique(clustering))
    clusters <- setdiff(clusters, 0) #Remove the noise cluster
    
    node <- lapply(clusters, function(cluster) {
      is_cluster <- clustering == cluster
      size <- sum(is_cluster)
      
      res_clean <- as.numeric(gsub(prefix, "", res))
      node_name <- paste0(prefix, res_clean, "C", cluster)
      
      node_data <- list(node_name, res_clean, cluster, size)
      names(node_data) <- c("node", prefix, "cluster", "size")
      
      for (aes in node_aes_list) {
        node_data <- clustree:::aggr_metadata(node_data, aes[[1]], aes[[2]],
                                   metadata, is_cluster)
      }
      
      node_data <- data.frame(node_data, stringsAsFactors = FALSE)
      
      return(node_data)
    })
    
    
    node <- do.call("rbind", node)
    
  })
  
  nodes <- do.call("rbind", nodes)
  
  return(nodes)
}


## Code modified from the
## clustertree (https://github.com/lazappi/clustree) package:

#' Get tree edges
#'
#' Extract the edges from a set of clusterings
#'
#' @param clusterings numeric matrix containing clustering information, each
#' column contains clustering at a separate resolution
#' @param prefix string indicating columns containing clustering information
#'
#' @return data.frame containing edge information
#'
#' @keywords internal
#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data
get_tree_edges_p <- function(clusterings, prefix) {
  
  res_values <- colnames(clusterings)
  
  edges <- lapply(seq_len(ncol(clusterings) - 1), function(idx) {
    from_res <- res_values[idx]
    to_res <- res_values[idx + 1]
    
    #Remove the zero cluster
    from_clusters <- setdiff(sort(unique(clusterings[, from_res])), 0)
    to_clusters <- setdiff(sort(unique(clusterings[, to_res])), 0)
    
    from_tos <- expand.grid(from_clust = from_clusters,
                            to_clust = to_clusters,
                            stringsAsFactors = FALSE)
    
    transitions <- apply(from_tos, 1, function(from_to) {
      from_clust <- from_to[1]
      to_clust <- from_to[2]
      
      is_from <- clusterings[, from_res] == from_clust
      is_to <- clusterings[, to_res] == to_clust
      
      trans_count <- sum(is_from & is_to)
      
      to_size <- sum(is_to)
      
      in_prop <- trans_count / to_size
      
      return(c(trans_count, in_prop))
    })
    
    from_tos$from_res <- as.numeric(gsub(prefix, "", from_res))
    from_tos$to_res <- as.numeric(gsub(prefix, "", to_res))
    from_tos$count <- transitions[1, ]
    from_tos$in_prop <- transitions[2, ]
    
    return(from_tos)
  })
  
  edges <- dplyr::bind_rows(edges) %>%
    dplyr::mutate(from_node = paste0(prefix, !!as.name("from_res"),
                                     "C", !!as.name("from_clust"))) %>%
    dplyr::mutate(to_node = paste0(prefix, !!as.name("to_res"),
                                   "C", !!as.name("to_clust"))) %>%
    dplyr::select("from_node", "to_node", dplyr::everything()) %>%
    dplyr::rename(!!as.name(paste0("from_", prefix)) := !!as.name("from_res"),
                  !!as.name(paste0("to_", prefix)) := !!as.name("to_res"))
  
  return(edges)
  
}

if(any(grepl("package:clustree", search()))) {
  message("Clustree deteched")
  detach("package:clustree") 
}

library(clustree)

assignInNamespace("get_tree_nodes",
                  get_tree_nodes_p,
                  ns = "clustree")
assignInNamespace("get_tree_edges",
                  get_tree_edges_p,
                  ns = "clustree")