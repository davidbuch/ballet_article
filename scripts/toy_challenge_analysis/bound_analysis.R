library(ggplot2)
library(grid)
library(gridExtra)
library(dirichletprocess)
source("R/dirichletprocess_custom_init.R")
source("R/rearrange_labels.R")
Rcpp::sourceCpp("src/subpartiton_min_max.cpp")
dir.create('output/toy_challenge', recursive = TRUE, showWarnings = FALSE)
set.seed(12345)

# Helper functions
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
nunique <- function(lbs){
  length(unique(lbs))
}
color_vals_ts <- c(gg_color_hue(10), 
                   rep(RColorBrewer::brewer.pal(4, 'Greens'), 1000))


# Load Data and Set Parameters
x <- scale(as.matrix(readRDS("data/clean_data/tsne.rds")))
nobs <- nrow(x)
nsims <- 1000
alpha <- 0.05

# Fit the DPMM  Model
g0Priors <- list(mu0 = rep(0,ncol(x)),
                 Lambda = diag(ncol(x)),
                 kappa0 = 0.1, #0.01,
                 nu = 10) # 200
dp_mod <- DirichletProcessMvnormal(x, 
                                   numInitialClusters = 20,
                                   g0Priors = g0Priors,
                                   alphaPriors = c(1,0.01))
dp_mod <- custom_init(dp_mod)
dp_mod <- Fit(dp_mod, nsims)

# extract mcmc samples of f
mixture_clustering_samps <- matrix(nrow = nsims, ncol = nobs)
for(s in 1:nsims)
  mixture_clustering_samps[s,] <- dp_mod$labelsChain[[s]]

# Find the Clustering Point Estimate (VI Loss 1,1)
mixture_pe <- salso::salso(mixture_clustering_samps,
                           loss = salso::VI(),
                           maxNClusters = 1000,
                           maxZealousAttempts = 1000)
nunique(mixture_pe)
table(apply(mixture_clustering_samps, 1, nunique))

# Precompute relevant quantities for point estimation and uncertainty bounds
psm <- mcclust::comp.psm(mixture_clustering_samps)
sample_losses <- apply(mixture_clustering_samps, 1, 
                       \(cs) salso::VI(mixture_pe, cs))
eps_star <- quantile(sample_losses, 1 - alpha)

# Method 1 - Wade and Gharamani Extreme Samples from Credible Ball
wg_bounds <- mcclust.ext::credibleball(mixture_pe, 
                                       mixture_clustering_samps, 
                                       c.dist = 'VI')
wglb <- wg_bounds$c.lowervert[1,]
wgub <- wg_bounds$c.uppervert[1,]

# Method 2 - Lattice Min and Lattice Max of Samples from Credible Ball
credible_ball <- mixture_clustering_samps[sample_losses < eps_star,]
min_of_ball <- min_subpartition(credible_ball)
max_of_ball <- max_subpartition(credible_ball)

# Method 3 - Find Low Risk partitions that fall just outside the Ball
# (Another Implementation would be to scan over salso and minimize risk conditional on fixed numbers of Clusters - I think the results would be similar)
# Scan over the value of Sep Loss
bp_grid <- seq(0.1, 1.9, by = 0.1)
mixture_pes <- matrix(nrow = length(bp_grid), ncol = ncol(mixture_clustering_samps))
for(pv in 1:nrow(mixture_pes)){
  cat(sprintf("Penalty Value: %d\n", pv))
  print(bench::bench_time(mixture_pes[pv,] <- salso::salso(mixture_clustering_samps, 
                                   loss = salso::VI(a = bp_grid[pv]),
                                   maxNClusters = 1000,
                                   maxZealousAttempts = 1000)))
}
pe_losses <- apply(mixture_pes, 1, \(cs) salso::VI(mixture_pe, cs))

png("output/toy_challenge/bound_analysis_greedy_distance.png", 
    width = 15, height = 20, units = 'in', res = 300)
plot(bp_grid - 1, pe_losses,
     ylim = c(max(0, min(pe_losses)), max(pe_losses, eps_star)),
     type = 'l', xlab = 'Bias Parameter', ylab = 'VI(1,1)', 
     main = 'Distance of Biased Partition\n from Point Estimate Partition')
abline(h = eps_star, col = 'red')
dev.off()

greedy_lb <- mixture_pes[min(which(pe_losses < eps_star)) - 1,]
greedy_ub <- mixture_pes[max(which(pe_losses < eps_star)) + 1,]


# Method 4 - Edgewise lower bound. Probably not coherent (the transitive completion of the edgewise lower bound might not be a lower bound of anythihng)
ewlb <- decode_subpartition(floyd_warshall(psm > 1 - alpha, rep(1,nobs)))
ewub <- decode_subpartition(floyd_warshall(psm > alpha, rep(1,nobs)))


# Rearrange labels to prepare for plotting
wglb <- rearrange_labels(wglb)
wgub <- rearrange_labels(wgub)

max_of_ball <- rearrange_labels(max_of_ball)
min_of_ball <- rearrange_labels(min_of_ball)

greedy_lb <- rearrange_labels(greedy_lb)
greedy_ub <- rearrange_labels(greedy_ub)

ewlb <- rearrange_labels(ewlb)
ewub <- rearrange_labels(ewub)

p1 <- ggplot() + 
  geom_point(aes(x = x[,1], y = x[,2], color = factor(mixture_pe)), size = 0.1) + 
  guides(color = 'none') + scale_color_manual(values = color_vals_ts) +
  labs(title = sprintf("Point Estimate (%d)", nunique(mixture_pe)), 
       x = NULL, y = NULL)

p2 <- ggplot() + 
  geom_point(aes(x = x[,1], y = x[,2], color = factor(wglb)), size = 0.1) + 
  guides(color = 'none') + scale_color_manual(values = color_vals_ts) +
  labs(title = sprintf("WG Vertical Lower Bound (%d)", nunique(wglb)), 
       x = NULL, y = NULL)
p3 <- ggplot() + 
  geom_point(aes(x = x[,1], y = x[,2], color = factor(wgub)), size = 0.1) + 
  guides(color = 'none') + scale_color_manual(values = color_vals_ts) + 
  labs(title = sprintf("WG Vertical Upper Bound (%d)", nunique(wgub)), 
       x = NULL, y = NULL)


p4 <- ggplot() + 
  geom_point(aes(x = x[,1], y = x[,2], color = factor(min_of_ball)), size = 0.1) + 
  guides(color = 'none') + scale_color_manual(values = color_vals_ts) +
  labs(title = sprintf("Lattice Min of Credible Ball (%d)", nunique(min_of_ball)), 
       x = NULL, y = NULL)
p5 <- ggplot() + 
  geom_point(aes(x = x[,1], y = x[,2], color = factor(max_of_ball)), size = 0.1) + 
  guides(color = 'none') + scale_color_manual(values = color_vals_ts) + 
  labs(title = sprintf("Lattice Max of Credible Ball (%d)", nunique(max_of_ball)), 
       x = NULL, y = NULL)

p6 <- ggplot() + 
  geom_point(aes(x = x[,1], y = x[,2], color = factor(greedy_lb)), size = 0.1) + 
  guides(color = 'none') + scale_color_manual(values = color_vals_ts) +
  labs(title = sprintf("Greedy Lower Bound of Credible Ball (%d)", 
                       nunique(greedy_lb)), 
       x = NULL, y = NULL)
p7 <- ggplot() + 
  geom_point(aes(x = x[,1], y = x[,2], color = factor(greedy_ub)), size = 0.1) + 
  guides(color = 'none') + scale_color_manual(values = color_vals_ts) + 
  labs(title = sprintf("Greedy Upper Bound of Credible Ball (%d)", 
                       nunique(greedy_ub)), 
       x = NULL, y = NULL)

p8 <- ggplot() + 
  geom_point(aes(x = x[,1], y = x[,2], color = factor(ewlb)), size = 0.1) + 
  guides(color = 'none') + scale_color_manual(values = color_vals_ts) + 
  labs(title = sprintf("Edge-wise Lower Bound (%d)", nunique(ewlb)), 
       x = NULL, y = NULL)
p9 <- ggplot() + 
  geom_point(aes(x = x[,1], y = x[,2], color = factor(ewub)), size = 0.1) + 
  guides(color = 'none') + scale_color_manual(values = color_vals_ts) + 
  labs(title = sprintf("Edge-wise Upper Bound (%d)", nunique(ewub)), 
       x = NULL, y = NULL)

png("output/toy_challenge/bound_analysis.png", 
    width = 15, height = 20, units = 'in', res = 300)
grid.arrange(p2, p1, p3,
             p4, p1, p5,
             p6, p1, p7, 
             p8, p1, p9,
             nrow = 4, ncol = 3, 
             top = textGrob("Bound Analysis - Model Based Clustering", 
                            gp = gpar(fontsize = 18)))
dev.off()
