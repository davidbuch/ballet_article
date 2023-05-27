# Fit the usual DBSCAN algorithm to the datasets
library(dbscan)

# minPts = 4 since dim = 2
dbscan::kNNdistplot(tsne,4)
abline(h = 2)

# fit with eps = 2 based on the k-dist plot
dbf <- dbscan::dbscan(tsne, 2, 4)

# remove non-prominent clusters
ccounts <- table(dbf$cluster)
clusters <- sapply(dbf$cluster, 
                   \(x) ifelse(ccounts[as.character(x)] < 0.01*nrow(tsne), 0, x))

# plot the tSNE DBSCAN fit
par(mfrow = c(1,2))
plot(tsne)
plot(tsne, col = clusters, pch = clusters)

