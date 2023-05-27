f <- readRDS("fitted_models/simulated_data/aniso_dpmm.rds")
x <- readRDS("data/clean_data/simulated/aniso.rds")
ggplot() + 
  geom_point(aes(x = x[,1], y = x[,2], color = colMeans(f)))

ggplot() + 
  geom_point(aes(x = x[,1], y = x[,2]), color = "grey50") + 
  geom_point(aes(x = x[f > tau,1], y = x[f > tau,2], color = factor(cl)))

ggplot() + 
  geom_point(aes(x = x[,1], y = x[,2], color = factor(labels)))


# Nothing about this theory bound for delta inherently requires that delta
# decrease with N... The idea is that the tube may just be \eta off from the 
# true tube... I suppose you could decrease eta? How quickly can you decrease eta?
# I think I am probably better off using my empirical method and later determining 
# whether I have set \delta < \delta_\epsilon

library(RcppAlgos)
f <- readRDS("fitted_models/simulated_data/aniso_dpmm.rds")
x <- readRDS("data/clean_data/simulated/aniso.rds")

N <- 500
sel <- sample(1:nrow(x), N, replace=TRUE)
x <- x[sel,]
f <- f[,sel]


D <- ncol(x)

minPts <- floor(0.05 * N)
fmean <- colMeans(f)
fsorted <- sort(fmean)

eta <- min(fsorted[(minPts + 1):N] - fsorted[1:(N - minPts)]) / 16
eta <- fsorted[N] - fsorted[round(0.15*N)]

xy <- cbind(x, fmean)
ind <- comboGeneral(nrow(xy), 2)
d <- xy[ind[, 2], ] - xy[ind[, 1], ]
slope <- d[, D + 1] / sqrt(rowSums(d[, 1:D, drop = FALSE]^2))
max(slope, na.rm = TRUE)

holder_const <- 2*max(slope, na.rm = TRUE)
holder_const

delta <- 0.5 * (eta / holder_const)^(1/alpha)
N
delta


