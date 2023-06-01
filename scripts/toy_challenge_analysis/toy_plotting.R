library(tidyverse)
library(grid)
library(gridExtra)
library(concaveman)
source("R/rearrange_labels.R")

# --------------------------------
# PART I: DPMM PLOTS
# --------------------------------
# Load the Plotting Objects
plot_obs_tm <- readRDS("output/toy_challenge/plot_obs_dpmm_two_moons.rds")
plot_obs_nc <- readRDS("output/toy_challenge/plot_obs_dpmm_circles.rds")
plot_obs_ts <- readRDS("output/toy_challenge/plot_obs_dpmm_tsne.rds")

plot_grid_tm <- readRDS("output/toy_challenge/plot_grid_dpmm_two_moons.rds")
plot_grid_nc <- readRDS("output/toy_challenge/plot_grid_dpmm_circles.rds")
plot_grid_ts <- readRDS("output/toy_challenge/plot_grid_dpmm_tsne.rds")

# Rearrange labels to facilitate better color discrimination among the largest clusters
prep_labels <- function(plot_obs){
  label_columns <- setdiff(colnames(plot_obs), c('x', 'y'))
  for(lc in label_columns){
    plot_obs[[lc]] <- rearrange_labels(plot_obs[[lc]])
    plot_obs[[lc]] <- factor(plot_obs[[lc]], 
                             levels = 1:max(plot_obs[[lc]]))
  }
  return(plot_obs)
}
plot_obs_tm <- prep_labels(plot_obs_tm)
plot_obs_nc <- prep_labels(plot_obs_nc)
plot_obs_ts <- prep_labels(plot_obs_ts)

# Plot coloring palette creation function
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
color_vals_nc <- c(gg_color_hue(2), 
                   rep(RColorBrewer::brewer.pal(4, 'Greens'), 1000))
color_vals_tm <- c(gg_color_hue(6), 
                   rep(RColorBrewer::brewer.pal(4, 'Greens'), 1000))
color_vals_ts <- c(gg_color_hue(10), 
                   rep(RColorBrewer::brewer.pal(4, 'Greens'), 1000))

# Function to help us set plot titles based on cluster richness
ktitle <- function(x) {
  sprintf("K = %d (%d)", nlevels(x), sum(table(x) > 1))
}

# First we will visualize:
# - Density contours and heatplots
# - Model based cluster allocations
# - Density based cluster allocations
pf_tm <- ggplot(plot_grid_tm) +   
  geom_contour_filled(aes(x = x, y = y, z = f_pe)) + 
  guides(fill = 'none') + 
  labs(x = NULL, y = NULL)
pf_ts <- ggplot(plot_grid_ts) +   
  geom_contour_filled(aes(x = x, y = y, z = f_pe)) + 
  guides(fill = 'none') + 
  labs(x = NULL, y = NULL)

pmc_tm <- ggplot(plot_obs_tm) + 
  geom_point(aes(x = x, y = y, color = mm_pe), size = 0.5) + 
  stat_ellipse(aes(x = x, y = y, group = mm_pe)) + 
  scale_color_manual(values = color_vals_tm) +
  guides(color = 'none') +
  theme(plot.title = element_text(size=9)) +
  labs(title = ktitle(plot_obs_tm$mm_pe), x = NULL, y = NULL)
pmc_ts <- ggplot(plot_obs_ts) + 
  geom_point(aes(x = x, y = y, color = mm_pe), size = 0.25) + 
  stat_ellipse(aes(x = x, y = y, group = mm_pe)) + 
  scale_color_manual(values = color_vals_ts) +
  guides(color = 'none') + 
  theme(plot.title = element_text(size=9)) +
  labs(title = ktitle(plot_obs_ts$mm_pe), x = NULL, y = NULL)

pdc_tm <- ggplot() + 
  geom_point(aes(x = x, y = y, color = db_pe), 
             size = 0.5, data = plot_obs_tm) + 
  geom_contour(aes(x = x, y = y, z = f_pe), 
               breaks = 0.09, col = 'black',
               data = plot_grid_tm) + 
  scale_color_manual(values = color_vals_tm) +
  guides(color = 'none') + 
  theme(plot.title = element_text(size=9)) +
  labs(title = ktitle(plot_obs_tm$db_pe), x = NULL, y = NULL)
pdc_ts <- ggplot() + 
  geom_point(aes(x = x, y = y, color = db_pe), 
             size = 0.25, data = plot_obs_ts) + 
  ggforce::geom_mark_hull(
    aes(x = x, y = y, group = db_pe), 
    expand = 1e-2,
    radius = 1e-2,
    concavity = 2,
    data = plot_obs_ts %>% filter(db_pe != 0)) +
  scale_color_manual(values = color_vals_ts) +
  guides(color = 'none') + 
  theme(plot.title = element_text(size=9)) +
  labs(title = ktitle(plot_obs_ts$db_pe), x = NULL, y = NULL)

png("output/toy_challenge/compare_clusterings_dpmm.png", 
    width = 12, height = 8, units = 'in', res = 300)
grid.arrange(
  arrangeGrob(pf_tm, pf_ts, 
              top = textGrob(
                "DPMM Density Estimate",
                gp = gpar(fontface = 3, fontsize = 14)
              )),
  arrangeGrob(pmc_tm, pmc_ts, 
              top = textGrob(
                "Model-Based",
                gp = gpar(fontface = 3, fontsize = 14)
              )),
  arrangeGrob(pdc_tm, pdc_ts,
              top = textGrob(
                "BAND",
                gp = gpar(fontface = 3, fontsize = 14)
              )),
  top = textGrob("Clustering Point Estimates - DPMM", 
                 gp = gpar(fontsize = 18)),
  ncol = 3
)
dev.off()

# Then we will visualize model-based clusterings and 
# credible bounds for each dataset
pvl_tm <- ggplot(plot_obs_tm) + 
  geom_point(aes(x = x, y = y, color = mm_vl_wg), size = 0.5) + 
  stat_ellipse(aes(x = x, y = y, group = mm_vl_wg)) + 
  scale_color_manual(values = color_vals_tm) +
  guides(color = 'none') + 
  theme(plot.title = element_text(size=9)) +
  labs(title = ktitle(plot_obs_tm$mm_vl_wg), x = NULL, y = NULL)
pvl_ts <- ggplot(plot_obs_ts) + 
  geom_point(aes(x = x, y = y, color = mm_vl_wg), size = 0.25) + 
  stat_ellipse(aes(x = x, y = y, group = mm_vl_wg)) + 
  scale_color_manual(values = color_vals_ts) +
  guides(color = 'none') + 
  theme(plot.title = element_text(size=9)) +
  labs(title = ktitle(plot_obs_ts$mm_vl_wg), x = NULL, y = NULL)

pvu_tm <- ggplot(plot_obs_tm) + 
  geom_point(aes(x = x, y = y, color = mm_vu_wg), size = 0.5) + 
  stat_ellipse(aes(x = x, y = y, group = mm_vu_wg)) + 
  scale_color_manual(values = color_vals_tm) +
  guides(color = 'none') + 
  theme(plot.title = element_text(size=9)) +
  labs(title = ktitle(plot_obs_tm$mm_vu_wg), x = NULL, y = NULL)
pvu_ts <- ggplot(plot_obs_ts) + 
  geom_point(aes(x = x, y = y, color = mm_vu_wg), size = 0.25) + 
  stat_ellipse(aes(x = x, y = y, group = mm_vu_wg)) + 
  scale_color_manual(values = color_vals_ts) +
  guides(color = 'none') + 
  theme(plot.title = element_text(size=9)) +
  labs(title = ktitle(plot_obs_ts$mm_vu_wg), x = NULL, y = NULL)

png("output/toy_challenge/mbc_wg_bounds_dpmm.png", 
    width = 12, height = 8, units = 'in', res = 300)
grid.arrange(
  arrangeGrob(pvl_tm, pvl_ts, 
              top = textGrob(
                "Lower Bound",
                gp = gpar(fontface = 3, fontsize = 14)
              )),
  arrangeGrob(pmc_tm, pmc_ts, 
              top = textGrob(
                "Point Estimate",
                gp = gpar(fontface = 3, fontsize = 14)
              )),
  arrangeGrob(pvu_tm, pvu_ts,
              top = textGrob(
                "Upper Bound",
                gp = gpar(fontface = 3, fontsize = 14)
              )),
  top = textGrob("Model-Based Clustering - W&G Credible Balls", 
                 gp = gpar(fontsize = 18)),
  ncol = 3
)
dev.off()

# Then we will visualize model-based clusterings (NOT W&G VERSION) and 
# credible bounds for each dataset
pvl_tm <- ggplot(plot_obs_tm) + 
  geom_point(aes(x = x, y = y, color = mm_vl), size = 0.5) + 
  stat_ellipse(aes(x = x, y = y, group = mm_vl)) + 
  scale_color_manual(values = color_vals_tm) +
  guides(color = 'none') + 
  theme(plot.title = element_text(size=9)) +
  labs(title = ktitle(plot_obs_tm$mm_vl), x = NULL, y = NULL)
pvl_ts <- ggplot(plot_obs_ts) + 
  geom_point(aes(x = x, y = y, color = mm_vl), size = 0.25) + 
  stat_ellipse(aes(x = x, y = y, group = mm_vl)) + 
  scale_color_manual(values = color_vals_ts) +
  guides(color = 'none') + 
  theme(plot.title = element_text(size=9)) +
  labs(title = ktitle(plot_obs_ts$mm_vl), x = NULL, y = NULL)

pvu_tm <- ggplot(plot_obs_tm) + 
  geom_point(aes(x = x, y = y, color = mm_vu), size = 0.5) + 
  stat_ellipse(aes(x = x, y = y, group = mm_vu)) + 
  scale_color_manual(values = color_vals_tm) +
  guides(color = 'none') + 
  theme(plot.title = element_text(size=9)) +
  labs(title = ktitle(plot_obs_tm$mm_vu), x = NULL, y = NULL)
pvu_ts <- ggplot(plot_obs_ts) + 
  geom_point(aes(x = x, y = y, color = mm_vu), size = 0.25) + 
  stat_ellipse(aes(x = x, y = y, group = mm_vu)) + 
  scale_color_manual(values = color_vals_ts) +
  guides(color = 'none') + 
  theme(plot.title = element_text(size=9)) +
  labs(title = ktitle(plot_obs_ts$mm_vu), x = NULL, y = NULL)

png("output/toy_challenge/mbc_bounds_dpmm.png", 
    width = 12, height = 8, units = 'in', res = 300)
grid.arrange(
  arrangeGrob(pvl_tm, pvl_ts, 
              top = textGrob(
                "Lower Bound",
                gp = gpar(fontface = 3, fontsize = 14)
              )),
  arrangeGrob(pmc_tm, pmc_ts, 
              top = textGrob(
                "Point Estimate",
                gp = gpar(fontface = 3, fontsize = 14)
              )),
  arrangeGrob(pvu_tm, pvu_ts,
              top = textGrob(
                "Upper Bound",
                gp = gpar(fontface = 3, fontsize = 14)
              )),
  top = textGrob("Model-Based Clustering - New Credible Balls", 
                 gp = gpar(fontsize = 18)),
  ncol = 3
)
dev.off()

# Finally we will visualize density-based clusterings and 
# credible bounds for each dataset
pvl_tm <- ggplot() + 
  geom_point(aes(x = x, y = y, color = db_vl),
             size = 0.5, data = plot_obs_tm) + 
  geom_contour(aes(x = x, y = y, z = f_pe), 
               breaks = 0.09, col = 'black',
               data = plot_grid_tm) + 
  scale_color_manual(values = color_vals_tm) +
  guides(color = 'none') + 
  theme(plot.title = element_text(size=9)) +
  labs(title = ktitle(plot_obs_tm$db_vl), x = NULL, y = NULL)

pvl_ts <- ggplot() + 
  geom_point(aes(x = x, y = y, color = db_vl), 
             size = 0.25, data = plot_obs_ts) + 
  ggforce::geom_mark_hull(
    aes(x = x, y = y, group = db_vl), 
    expand = 1e-2,
    radius = 1e-2,
    concavity = 2,
    data = plot_obs_ts %>% filter(db_vl != 0)) +
  scale_color_manual(values = color_vals_ts) + 
  guides(color = 'none') + 
  theme(plot.title = element_text(size=9)) +
  labs(title = ktitle(plot_obs_ts$db_vl), x = NULL, y = NULL)

pvu_tm <- ggplot() + 
  geom_point(aes(x = x, y = y, color = db_vu),
             size = 0.5, data = plot_obs_tm) + 
  geom_contour(aes(x = x, y = y, z = f_pe), 
               breaks = 0.09, col = 'black',
               data = plot_grid_tm) +
  scale_color_manual(values = color_vals_tm) + 
  guides(color = 'none') + 
  theme(plot.title = element_text(size=9)) +
  labs(title = ktitle(plot_obs_tm$db_vu), x = NULL, y = NULL)
pvu_ts <- ggplot() + 
  geom_point(aes(x = x, y = y, color = db_vu), 
             size = 0.25, data = plot_obs_ts) + 
  ggforce::geom_mark_hull(
    aes(x = x, y = y, group = db_vu), 
    expand = 1e-2,
    radius = 1e-2,
    concavity = 2,
    data = plot_obs_ts %>% filter(db_vu != 0)) +
  scale_color_manual(values = color_vals_ts) + 
  guides(color = 'none') + 
  theme(plot.title = element_text(size=9)) +
  labs(title = ktitle(plot_obs_ts$db_vu), x = NULL, y = NULL)

png("output/toy_challenge/band_bounds_dpmm.png", 
    width = 12, height = 8, units = 'in', res = 300)
grid.arrange(
  arrangeGrob(pvl_tm, pvl_ts, 
              top = textGrob(
                "Lower Bound",
                gp = gpar(fontface = 3, fontsize = 14)
              )),
  arrangeGrob(pdc_tm, pdc_ts, 
              top = textGrob(
                "Point Estimate",
                gp = gpar(fontface = 3, fontsize = 14)
              )),
  arrangeGrob(pvu_tm, pvu_ts,
              top = textGrob(
                "Upper Bound",
                gp = gpar(fontface = 3, fontsize = 14)
              )),
  top = textGrob("BAND Clustering - Credible Balls", 
                 gp = gpar(fontsize = 18)),
  ncol = 3
)
dev.off()

# -----------------------------------
# PART II: COMPARING DENSITY MODELS
# -----------------------------------
plot_grid_tm_apt <- readRDS("output/toy_challenge/plot_grid_apt_two_moons.rds")
plot_grid_tm_nndm<- readRDS("output/toy_challenge/plot_grid_nndm_two_moons.rds")
plot_grid_ts_apt <- readRDS("output/toy_challenge/plot_grid_apt_tsne.rds")
plot_grid_ts_nndm <- readRDS("output/toy_challenge/plot_grid_nndm_tsne.rds")
plot_grid_nc_apt <- readRDS("output/toy_challenge/plot_grid_apt_circles.rds")
plot_grid_nc_nndm<- readRDS("output/toy_challenge/plot_grid_nndm_circles.rds")

plot_obs_tm_apt <- prep_labels(readRDS("output/toy_challenge/plot_obs_apt_two_moons.rds"))
plot_obs_tm_nndm<- prep_labels(readRDS("output/toy_challenge/plot_obs_nndm_two_moons.rds"))
plot_obs_ts_apt <- prep_labels(readRDS("output/toy_challenge/plot_obs_apt_tsne.rds"))
plot_obs_ts_nndm <- prep_labels(readRDS("output/toy_challenge/plot_obs_nndm_tsne.rds"))
plot_obs_nc_apt <- prep_labels(readRDS("output/toy_challenge/plot_obs_apt_circles.rds"))
plot_obs_nc_nndm<- prep_labels(readRDS("output/toy_challenge/plot_obs_nndm_circles.rds"))

# Density Estimate Plots (there should be 9 of them)
pf_tm_dpmm <- ggplot(plot_grid_tm) +   
  geom_contour_filled(aes(x = x, y = y, z = f_pe)) + 
  guides(fill = 'none') + 
  labs(x = NULL, y = NULL)

pf_tm_apt <- ggplot(plot_grid_tm_apt) +   
  geom_contour_filled(aes(x = x, y = y, z = f_pe)) + 
  guides(fill = 'none') + 
  labs(x = NULL, y = NULL)

pf_tm_nndm <- ggplot(plot_grid_tm_nndm) +   
  geom_contour_filled(aes(x = x, y = y, z = f_pe)) + 
  guides(fill = 'none') + 
  labs(x = NULL, y = NULL)

pf_ts_dpmm <- ggplot(plot_grid_ts) +   
  geom_contour_filled(aes(x = x, y = y, z = f_pe)) + 
  guides(fill = 'none') + 
  labs(x = NULL, y = NULL)

pf_ts_apt <- ggplot(plot_grid_ts_apt) +   
  geom_contour_filled(aes(x = x, y = y, z = f_pe)) + 
  guides(fill = 'none') + 
  labs(x = NULL, y = NULL)

pf_ts_nndm <- ggplot(plot_grid_ts_nndm) +   
  geom_contour_filled(aes(x = x, y = y, z = f_pe)) + 
  guides(fill = 'none') + 
  labs(x = NULL, y = NULL)

pf_nc_dpmm <- ggplot(plot_grid_nc) +   
  geom_contour_filled(aes(x = x, y = y, z = f_pe)) + 
  guides(fill = 'none') + 
  labs(x = NULL, y = NULL)

pf_nc_apt <- ggplot(plot_grid_nc_apt) +   
  geom_contour_filled(aes(x = x, y = y, z = f_pe)) + 
  guides(fill = 'none') + 
  labs(x = NULL, y = NULL)

pf_nc_nndm <- ggplot(plot_grid_nc_nndm) +   
  geom_contour_filled(aes(x = x, y = y, z = f_pe)) + 
  guides(fill = 'none') + 
  labs(x = NULL, y = NULL)

png("output/toy_challenge/compare_density_models.png", 
    width = 12, height = 8, units = 'in', res = 300)
grid.arrange(
  arrangeGrob(pf_tm_dpmm, pf_nc_dpmm, pf_ts_dpmm, 
              top = textGrob(
                "DP Mixture of Gaussians",
                gp = gpar(fontface = 3, fontsize = 14)
              )),
  arrangeGrob(pf_tm_apt, pf_nc_apt, pf_ts_apt,
              top = textGrob(
                "Adaptive Polya Tree",
                gp = gpar(fontface = 3, fontsize = 14)
              )),
  arrangeGrob(pf_tm_nndm, pf_nc_nndm, pf_ts_nndm,
              top = textGrob(
                "NN Dirichlet Mixture",
                gp = gpar(fontface = 3, fontsize = 14)
              )),
  top = textGrob("Density Model Point Estimates", 
                 gp = gpar(fontsize = 18)),
  ncol = 3
)
dev.off()

# Density clustering plots (there should be 9 of them)
pdc_tm_dpmm <- ggplot() + 
  geom_point(aes(x = x, y = y, color = db_pe), 
             size = 0.5, data = plot_obs_tm) + 
  geom_contour(aes(x = x, y = y, z = f_pe), 
               breaks = 0.09, col = 'black',
               data = plot_grid_tm) + 
  scale_color_manual(values = color_vals_tm) +
  guides(color = 'none') + 
  theme(plot.title = element_text(size=9)) +
  labs(title = ktitle(plot_obs_tm$db_pe), x = NULL, y = NULL)

pdc_nc_dpmm <- ggplot() + 
  geom_point(aes(x = x, y = y, color = db_pe), 
             size = 0.5, data = plot_obs_nc) + 
  geom_contour(aes(x = x, y = y, z = f_pe), 
               breaks = 0.1, col = 'black',
               data = plot_grid_nc) + 
  scale_color_manual(values = color_vals_nc) +
  guides(color = 'none') + 
  theme(plot.title = element_text(size=9)) +
  labs(title = ktitle(plot_obs_nc$db_pe), x = NULL, y = NULL)

pdc_ts_dpmm <- ggplot() + 
  geom_point(aes(x = x, y = y, color = db_pe), 
             size = 0.25, data = plot_obs_ts) + 
  ggforce::geom_mark_hull(
    aes(x = x, y = y, group = db_pe), 
    expand = 1e-2,
    radius = 1e-2,
    concavity = 2,
    data = plot_obs_ts %>% filter(db_pe != 0)) +
  scale_color_manual(values = color_vals_ts) +
  guides(color = 'none') + 
  theme(plot.title = element_text(size=9)) +
  labs(title = ktitle(plot_obs_ts$db_pe), x = NULL, y = NULL)

pdc_tm_apt <- ggplot() + 
  geom_point(aes(x = x, y = y, color = db_pe), 
             size = 0.5, data = plot_obs_tm_apt) + 
  geom_contour(aes(x = x, y = y, z = f_pe), 
               breaks = 0.09, col = 'black',
               data = plot_grid_tm_apt) + 
  scale_color_manual(values = color_vals_tm) +
  guides(color = 'none') + 
  theme(plot.title = element_text(size=9)) +
  labs(title = ktitle(plot_obs_tm_apt$db_pe), x = NULL, y = NULL)

pdc_nc_apt <- ggplot() + 
  geom_point(aes(x = x, y = y, color = db_pe), 
             size = 0.5, data = plot_obs_nc_apt) + 
  geom_contour(aes(x = x, y = y, z = f_pe), 
               breaks = 0.07, #quantile(plot_obs_nc_apt$f_pe, 0.125),
               col = 'black',
               data = plot_grid_nc_apt) + 
  scale_color_manual(values = color_vals_tm) +
  guides(color = 'none') + 
  theme(plot.title = element_text(size=9)) +
  labs(title = ktitle(plot_obs_nc_apt$db_pe), x = NULL, y = NULL)

pdc_ts_apt <- ggplot() + 
  geom_point(aes(x = x, y = y, color = db_pe), 
             size = 0.25, data = plot_obs_ts_apt) + 
  ggforce::geom_mark_hull(
    aes(x = x, y = y, group = db_pe), 
    expand = 1e-2,
    radius = 1e-2,
    concavity = 2,
    data = plot_obs_ts_apt %>% filter(db_pe != 0)) +
  scale_color_manual(values = color_vals_ts) +
  guides(color = 'none') + 
  theme(plot.title = element_text(size=9)) +
  labs(title = ktitle(plot_obs_ts_apt$db_pe), x = NULL, y = NULL)

pdc_tm_nndm <- ggplot() + 
  geom_point(aes(x = x, y = y, color = db_pe), 
             size = 0.5, data = plot_obs_tm_nndm) + 
  geom_contour(aes(x = x, y = y, z = f_pe), 
               breaks = 0.11, col = 'black',
               data = plot_grid_tm_nndm) + 
  scale_color_manual(values = color_vals_tm) +
  guides(color = 'none') + 
  theme(plot.title = element_text(size=9)) +
  labs(title = ktitle(plot_obs_tm_nndm$db_pe), x = NULL, y = NULL)

pdc_nc_nndm <- ggplot() + 
  geom_point(aes(x = x, y = y, color = db_pe), 
             size = 0.5, data = plot_obs_nc_nndm) + 
  geom_contour(aes(x = x, y = y, z = f_pe), 
               breaks = 0.14,
               col = 'black',
               data = plot_grid_nc_nndm) + 
  scale_color_manual(values = color_vals_tm) +
  guides(color = 'none') + 
  theme(plot.title = element_text(size=9)) +
  labs(title = ktitle(plot_obs_nc_nndm$db_pe), x = NULL, y = NULL)

pdc_ts_nndm <- ggplot() + 
  geom_point(aes(x = x, y = y, color = db_pe), 
             size = 0.25, data = plot_obs_ts_nndm) + 
  ggforce::geom_mark_hull(
    aes(x = x, y = y, group = db_pe), 
    expand = 1e-2,
    radius = 1e-2,
    concavity = 2,
    data = plot_obs_ts_nndm %>% filter(db_pe != 0)) +
  scale_color_manual(values = color_vals_ts) +
  guides(color = 'none') + 
  theme(plot.title = element_text(size=9)) +
  labs(title = ktitle(plot_obs_ts_nndm$db_pe), x = NULL, y = NULL)

png("output/toy_challenge/compare_band_clusterings.png", 
    width = 12, height = 8, units = 'in', res = 300)
grid.arrange(
  arrangeGrob(pdc_tm_dpmm, pdc_nc_dpmm, pdc_ts_dpmm, 
              top = textGrob(
                "DP Mixture of Gaussians",
                gp = gpar(fontface = 3, fontsize = 14)
              )),
  arrangeGrob(pdc_tm_apt, pdc_nc_apt, pdc_ts_apt,
              top = textGrob(
                "Adaptive Polya Tree",
                gp = gpar(fontface = 3, fontsize = 14)
              )),
  arrangeGrob(pdc_tm_nndm, pdc_nc_nndm, pdc_ts_nndm,
              top = textGrob(
                "NN Dirichlet Mixture",
                gp = gpar(fontface = 3, fontsize = 14)
              )),
  top = textGrob("BAND Clustering Point Estimates", 
                 gp = gpar(fontsize = 18)),
  ncol = 3
)
dev.off()

# pf_nc <- ggplot(plot_grid_nc) +   
#   geom_contour_filled(aes(x = x, y = y, z = f_pe)) + 
#   guides(fill = 'none') + 
#   labs(x = NULL, y = NULL)
# pmc_nc <- ggplot(plot_obs_nc) + 
#   geom_point(aes(x = x, y = y, color = mm_pe)) + 
#   stat_ellipse(aes(x = x, y = y, group = mm_pe)) + 
#   guides(color = 'none') + 
#   labs(x = NULL, y = NULL)
# pdc_nc <- ggplot() + 
#   geom_point(aes(x = x, y = y, color = db_pe), 
#              data = plot_obs_nc) + 
#   geom_contour(aes(x = x, y = y, z = f_pe), 
#                breaks = 0.07, col = 'black',
#                data = plot_grid_nc) + 
#   guides(color = 'none') + 
#   labs(x = NULL, y = NULL)
# pdc_nc <- ggplot() + 
#   geom_point(aes(x = x, y = y, color = db_pe), data = plot_obs_nc) + 
#   ggforce::geom_mark_hull(
#     aes(x = x, y = y, group = db_pe), 
#     expand = 2e-2,
#     radius = 1e-2,
#     concavity = 2,
#     data = plot_obs_nc %>% filter(db_pe != 0)) +
#   guides(color = 'none') + 
#   labs(x = NULL, y = NULL)
# 
# pvl_nc <- ggplot() + 
#   geom_point(aes(x = x, y = y, color = db_vl), 
#              data = plot_obs_nc) + 
#   geom_contour(aes(x = x, y = y, z = f_pe), 
#                breaks = 0.07, col = 'black',
#                data = plot_grid_nc) + 
#   guides(color = 'none') + 
#   labs(x = NULL, y = NULL)
# pvu_nc <- ggplot() + 
#   geom_point(aes(x = x, y = y, color = db_vu), 
#              data = plot_obs_nc) + 
#   geom_contour(aes(x = x, y = y, z = f_pe), 
#                breaks = 0.07, col = 'black',
#                data = plot_grid_nc) + 
#   guides(color = 'none') + 
#   labs(x = NULL, y = NULL)
