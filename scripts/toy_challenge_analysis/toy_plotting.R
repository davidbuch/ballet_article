library(tidyverse)
library(grid)
library(gridExtra)
library(concaveman)
library(latex2exp)
source("R/rearrange_labels.R")
source("R/choose_lambda.R")

# --------------------------------
# PART I: DPMM PLOTS
# --------------------------------

split_df_by_levels <- function(plot_obs) {
  
  lvls <- str_match(colnames(plot_obs), 'db_pe_(.+)$')[,2]
  lvls <- lvls[!is.na(lvls)]
  
  # Analyze levels to decide on their ordering
  lvls_num <- as.numeric(lvls)
  # preset levels are in multiples of 0.05,
  # while an elbow level will not typically be..
  elbow_level <- which(floor(lvls_num*1000) %% 50 != 0)
  
  if(is.na(elbow_level)) {
    warning("Could not find elbow level in [", paste(lvls, collapse = ","),"]")
  } else {
    lvls_num[elbow_level] <- NA
  }
  
  lvls <- lvls[rank(-lvls_num, na.last = TRUE)]

  plot_obs_list <- rep(list(NULL), length(lvls))

  for(i in seq_along(lvls)) {
    suffix <- lvls[i]
    plot_obs_list[[i]] <- rename_with(
      plot_obs,
      .fn = \(colnames) str_remove(colnames, paste0("_",suffix)),
      .cols = ends_with(suffix))
  }
  
  #names(plot_obs_list) <- lvls
  
  plot_obs_list
}

# Rearrange labels to facilitate better color discrimination among the largest clusters
prep_labels <- function(plot_obs) {
  label_columns <- setdiff(colnames(plot_obs), c('x', 'y', 'f_pe'))
  for(lc in label_columns){
    
    noise_frac <- attr(plot_obs[[lc]], 'noise_frac')
    
    plot_obs[[lc]] <- rearrange_labels(plot_obs[[lc]])
    plot_obs[[lc]] <- factor(plot_obs[[lc]], 
                             levels = 1:max(plot_obs[[lc]]))
    
    if(!is.null(noise_frac)) {
      attr(plot_obs[[lc]], 'noise_frac') <- noise_frac
    }
  }
  return(plot_obs)
}

# The final pre-processing
process <- function(plot_obs) {
  
  # separate the elbow and non-elbow estimates.
  plot_obs |> 
    split_df_by_levels() |> 
      map(prep_labels) 
}

# Load the Plotting Objects

readRDS("output/toy_challenge/plot_obs_dpmm_two_moons.rds") |>
  process() -> plot_obs_tm_ls
plot_obs_tm <- plot_obs_tm_ls[[1]]

readRDS("output/toy_challenge/plot_obs_dpmm_circles.rds") |>
  process() -> plot_obs_nc_ls
plot_obs_nc <- plot_obs_nc_ls[[1]]

readRDS("output/toy_challenge/plot_obs_dpmm_tsne.rds") |>
  process() -> plot_obs_ts_ls
plot_obs_ts <- plot_obs_ts_ls[[1]]


plot_grid_tm <- readRDS("output/toy_challenge/plot_grid_dpmm_two_moons.rds")
plot_grid_nc <- readRDS("output/toy_challenge/plot_grid_dpmm_circles.rds")
plot_grid_ts <- readRDS("output/toy_challenge/plot_grid_dpmm_tsne.rds")


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
# and noise fraction (if provided)
ktitle <- function(x) {
  if (is.null(attr(x,"noise_frac"))) {
    extra = ""
  } else {
    if (is.na(attr(x,"noise_frac"))) {
      extra = "Persistent clusters"
    } else {
      extra = sprintf("noise level: %.1f%%",100*attr(x,"noise_frac"))
    }
  }
  sprintf("K = %d (%d) %s", nlevels(x), sum(table(x) > 1), extra)
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
  scale_color_manual(na.value = "black", values = color_vals_tm) +
  guides(color = 'none') +
  theme(plot.title = element_text(size=9)) +
  labs(title = ktitle(plot_obs_tm$mm_pe), x = NULL, y = NULL)

pmc_ts <- ggplot(plot_obs_ts) + 
  geom_point(aes(x = x, y = y, color = mm_pe), size = 0.25) + 
  stat_ellipse(aes(x = x, y = y, group = mm_pe)) + 
  scale_color_manual(na.value = "black", values = color_vals_ts) +
  guides(color = 'none') + 
  theme(plot.title = element_text(size=9)) +
  labs(title = ktitle(plot_obs_ts$mm_pe), x = NULL, y = NULL)

pdc_tm <- ggplot() + 
  geom_point(aes(x = x, y = y, color = db_pe), 
             size = 0.5, data = plot_obs_tm) + 
  geom_contour(aes(x = x, y = y, z = f_pe), 
               breaks = 0.09, col = 'black',
               data = plot_grid_tm) + 
  scale_color_manual(na.value = "black", values = color_vals_tm) +
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
  scale_color_manual(na.value = "black", values = color_vals_ts) +
  guides(color = 'none') + 
  theme(plot.title = element_text(size=9)) +
  labs(title = ktitle(plot_obs_ts$db_pe), 
       x = NULL, y = NULL)

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
                "BALLET",
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
  scale_color_manual(na.value = "black", values = color_vals_tm) +
  guides(color = 'none') + 
  theme(plot.title = element_text(size=9)) +
  labs(title = ktitle(plot_obs_tm$mm_vl_wg), x = NULL, y = NULL)
pvl_ts <- ggplot(plot_obs_ts) + 
  geom_point(aes(x = x, y = y, color = mm_vl_wg), size = 0.25) + 
  stat_ellipse(aes(x = x, y = y, group = mm_vl_wg)) + 
  scale_color_manual(na.value = "black", values = color_vals_ts) +
  guides(color = 'none') + 
  theme(plot.title = element_text(size=9)) +
  labs(title = ktitle(plot_obs_ts$mm_vl_wg), x = NULL, y = NULL)

pvu_tm <- ggplot(plot_obs_tm) + 
  geom_point(aes(x = x, y = y, color = mm_vu_wg), size = 0.5) + 
  stat_ellipse(aes(x = x, y = y, group = mm_vu_wg)) + 
  scale_color_manual(na.value = "black", values = color_vals_tm) +
  guides(color = 'none') + 
  theme(plot.title = element_text(size=9)) +
  labs(title = ktitle(plot_obs_tm$mm_vu_wg), x = NULL, y = NULL)
pvu_ts <- ggplot(plot_obs_ts) + 
  geom_point(aes(x = x, y = y, color = mm_vu_wg), size = 0.25) + 
  stat_ellipse(aes(x = x, y = y, group = mm_vu_wg)) + 
  scale_color_manual(na.value = "black", values = color_vals_ts) +
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
  top = textGrob("Model-Based Clustering - W&G Credible Bounds", 
                 gp = gpar(fontsize = 18)),
  ncol = 3
)
dev.off()

# Then we will visualize model-based clusterings (NOT W&G VERSION) and 
# credible bounds for each dataset
pvl_tm <- ggplot(plot_obs_tm) + 
  geom_point(aes(x = x, y = y, color = mm_vl), size = 0.5) + 
  stat_ellipse(aes(x = x, y = y, group = mm_vl)) + 
  scale_color_manual(na.value = "black", values = color_vals_tm) +
  guides(color = 'none') + 
  theme(plot.title = element_text(size=9)) +
  labs(title = ktitle(plot_obs_tm$mm_vl), x = NULL, y = NULL)
pvl_ts <- ggplot(plot_obs_ts) + 
  geom_point(aes(x = x, y = y, color = mm_vl), size = 0.25) + 
  stat_ellipse(aes(x = x, y = y, group = mm_vl)) + 
  scale_color_manual(na.value = "black", values = color_vals_ts) +
  guides(color = 'none') + 
  theme(plot.title = element_text(size=9)) +
  labs(title = ktitle(plot_obs_ts$mm_vl), x = NULL, y = NULL)

pvu_tm <- ggplot(plot_obs_tm) + 
  geom_point(aes(x = x, y = y, color = mm_vu), size = 0.5) + 
  stat_ellipse(aes(x = x, y = y, group = mm_vu)) + 
  scale_color_manual(na.value = "black", values = color_vals_tm) +
  guides(color = 'none') + 
  theme(plot.title = element_text(size=9)) +
  labs(title = ktitle(plot_obs_tm$mm_vu), x = NULL, y = NULL)
pvu_ts <- ggplot(plot_obs_ts) + 
  geom_point(aes(x = x, y = y, color = mm_vu), size = 0.25) + 
  stat_ellipse(aes(x = x, y = y, group = mm_vu)) + 
  scale_color_manual(na.value = "black", values = color_vals_ts) +
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
  top = textGrob("Model-Based Clustering - New Credible Bounds", 
                 gp = gpar(fontsize = 18)),
  ncol = 3
)
dev.off()

# Finally we will visualize density-based clusterings and 
# credible bounds for each dataset

plot_obj_ballet_bounds_two_moons <- function(lvl=1) {
  
  plot_obs_tm <- plot_obs_tm_ls[[lvl]]
  
  pvl_tm <- ggplot() + 
    geom_point(aes(x = x, y = y, color = db_vl),
               size = 0.5, data = plot_obs_tm) + 
    geom_contour(aes(x = x, y = y, z = f_pe), 
                 breaks = 0.09, col = 'black',
                 data = plot_grid_tm) + 
    scale_color_manual(na.value = "black", values = color_vals_tm) +
    guides(color = 'none') + 
    theme(plot.title = element_text(size=9)) +
    labs(title = ktitle(plot_obs_tm$db_vl), x = NULL, y = NULL)


  pvu_tm <- ggplot() + 
    geom_point(aes(x = x, y = y, color = db_vu),
               size = 0.5, data = plot_obs_tm) + 
    geom_contour(aes(x = x, y = y, z = f_pe), 
                 breaks = 0.09, col = 'black',
                 data = plot_grid_tm) +
    scale_color_manual(na.value = "black", values = color_vals_tm) + 
    guides(color = 'none') + 
    theme(plot.title = element_text(size=9)) +
    labs(title = ktitle(plot_obs_tm$db_vu), x = NULL, y = NULL)
  
  list(lower=pvl_tm, upper=pvu_tm)
}


p1 <- plot_obj_ballet_bounds_two_moons(1)
p2 <- plot_obj_ballet_bounds_two_moons(2)

png(sprintf("output/toy_challenge/ballet_bounds_dpmm_two_moons.png"),
    width = 12, height = 8, units = 'in', res = 300)
grid.arrange(
  arrangeGrob(p1$lower, 
              top = textGrob(
                "Lower Bound",
                gp = gpar(fontface = 3, fontsize = 14)
              )),
  arrangeGrob(p1$upper,
              top = textGrob(
                "Upper Bound",
                gp = gpar(fontface = 3, fontsize = 14)
              )),
  arrangeGrob(p2$lower),
  arrangeGrob(p2$upper),
  top = textGrob("BALLET Clustering - Credible Bounds - Two Moons", 
                 gp = gpar(fontsize = 18)),
  ncol = 2
)
dev.off()


plot_obj_ballet_bounds_tsne <- function(lvl=1, file_suffix=sprintf("-%d",lvl)) {
  
  plot_obs_ts <- plot_obs_ts_ls[[lvl]]
  
  pvl_ts <- ggplot() + 
    geom_point(aes(x = x, y = y, color = db_vl), 
               size = 0.25, data = plot_obs_ts) + 
    ggforce::geom_mark_hull(
      aes(x = x, y = y, group = db_vl), 
      expand = 1e-2,
      radius = 1e-2,
      concavity = 2,
      data = plot_obs_ts %>% filter(db_vl != 0)) +
    scale_color_manual(na.value = "black", values = color_vals_ts) + 
    guides(color = 'none') + 
    theme(plot.title = element_text(size=9)) +
    labs(title = ktitle(plot_obs_ts$db_vl), x = NULL, y = NULL)
  
  pvu_ts <- ggplot() + 
    geom_point(aes(x = x, y = y, color = db_vu), 
               size = 0.25, data = plot_obs_ts) + 
    ggforce::geom_mark_hull(
      aes(x = x, y = y, group = db_vu), 
      expand = 1e-2,
      radius = 1e-2,
      concavity = 2,
      data = plot_obs_ts %>% filter(db_vu != 0)) +
    scale_color_manual(na.value = "black", values = color_vals_ts) + 
    guides(color = 'none') + 
    theme(plot.title = element_text(size=9)) +
    labs(title = ktitle(plot_obs_ts$db_vu), x = NULL, y = NULL)
  
  list(lower=pvl_ts, upper=pvu_ts)
}

p1 <- plot_obj_ballet_bounds_tsne(1)
p2 <- plot_obj_ballet_bounds_tsne(2)
p3 <- plot_obj_ballet_bounds_tsne(3)


png("output/toy_challenge/ballet_bounds_dpmm_tsne.png",
    width = 12, height = 4, units = 'in', res = 300)
grid.arrange(
  arrangeGrob(p1$lower, 
              top = textGrob(
                "Lower Bound",
                gp = gpar(fontface = 3, fontsize = 14)
              )),
  arrangeGrob(p1$upper,
              top = textGrob(
                "Upper Bound",
                gp = gpar(fontface = 3, fontsize = 14)
              )),
  top = textGrob("BALLET Clustering - Credible Bounds", 
                 gp = gpar(fontsize = 18)),
  ncol = 2
)
dev.off()

png("output/toy_challenge/ballet_bounds_dpmm_tsne_lvls.png",
    width = 12, height = 8, units = 'in', res = 300)
grid.arrange(
  arrangeGrob(p2$lower, 
              top = textGrob(
                "Lower Bound",
                gp = gpar(fontface = 3, fontsize = 14)
              )),
  arrangeGrob(p2$upper,
              top = textGrob(
                "Upper Bound",
                gp = gpar(fontface = 3, fontsize = 14)
              )),
  arrangeGrob(p3$lower),
  arrangeGrob(p3$upper),
  top = textGrob("BALLET Clustering - Credible Bounds", 
                 gp = gpar(fontsize = 18)),
  ncol = 2
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

readRDS("output/toy_challenge/plot_obs_apt_two_moons.rds") |> 
    process() -> plot_obs_tm_apt_ls

readRDS("output/toy_challenge/plot_obs_nndm_two_moons.rds") |>
    process() -> plot_obs_tm_nndm_ls

readRDS("output/toy_challenge/plot_obs_apt_tsne.rds") |>
    process() -> plot_obs_ts_apt_ls

readRDS("output/toy_challenge/plot_obs_nndm_tsne.rds") |>
    process() -> plot_obs_ts_nndm_ls

readRDS("output/toy_challenge/plot_obs_apt_circles.rds") |>
    process() -> plot_obs_nc_apt_ls

readRDS("output/toy_challenge/plot_obs_nndm_circles.rds") |>
    process() -> plot_obs_nc_nndm_ls

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

plot_func <- function(plot_obs, plot_grid, color_vals, 
                      size=0.5, break_thresh=0.125) {
  ggplot() + 
    geom_point(aes(x = x, y = y, color = db_pe), 
               size = size, data = plot_obs) + 
    geom_contour(aes(x = x, y = y, z = f_pe), 
                 breaks = quantile(plot_obs$f_pe, break_thresh), 
                 col = 'black',
                 data = plot_grid) + 
    scale_color_manual(na.value = "black", values = color_vals) +
    guides(color = 'none') + 
    theme(plot.title = element_text(size=9)) +
    labs(title = ktitle(plot_obs$db_pe), x = NULL, y = NULL)
}

plot_func_tsne <- function(plot_obs, size=0.25) {
    ggplot() + 
      geom_point(aes(x = x, y = y, color = db_pe), 
                 size = size, data = plot_obs) + 
      ggforce::geom_mark_hull(
        aes(x = x, y = y, group = db_pe), 
        expand = 1e-2,
        radius = 1e-2,
        concavity = 2,
        data = plot_obs %>% filter(db_pe != 0)) +
      scale_color_manual(na.value = "black", values = color_vals_ts) +
      guides(color = 'none') + 
      theme(plot.title = element_text(size=9)) +
      labs(title = ktitle(plot_obs$db_pe), x = NULL, y = NULL)
}

plot_across_density_models <- function(two_moons_lvl=1, circles_lvl=1, tsne_lvl=1,
                                      file_suffix=sprintf("%d-%d-%d", 
                                                      two_moons_lvl, 
                                                      circles_lvl,
                                                      tsne_lvl)) {
  
  safely_get_index <- function(ls, index) {
    N <- length(ls)
    
    if( index < 0 || index > N ) {
      warning("Index out of bounds. Safely truncating it.")
      index <- max(min(N, index), 0)
    }
    
    ls[[index]]
  }
  
  pdc_tm_dpmm <- plot_func(safely_get_index(plot_obs_tm_ls, two_moons_lvl), 
                           plot_grid_tm, color_vals_tm, size=0.5,  break_thresh=0.125)
  
  pdc_nc_dpmm <- plot_func(safely_get_index(plot_obs_nc_ls, circles_lvl), 
                           plot_grid_nc, color_vals_nc, size=0.5,  break_thresh=0.02)
  
  pdc_ts_dpmm <- plot_func_tsne(safely_get_index(plot_obs_ts_ls, tsne_lvl), size=0.25)
  
  pdc_tm_apt <- plot_func(safely_get_index(plot_obs_tm_apt_ls,two_moons_lvl), 
                          plot_grid_tm_apt, color_vals_tm,
                          size = 0.5, break_thresh = 0.08)
  
  pdc_nc_apt <- plot_func(safely_get_index(plot_obs_nc_apt_ls, circles_lvl), 
                          plot_grid_nc_apt, color_vals_nc, size = 0.5, 
                          break_thresh = 0.02)
  
  pdc_ts_apt <- plot_func_tsne(safely_get_index(plot_obs_ts_apt_ls, tsne_lvl), size=0.25)
  
  
  pdc_tm_nndm <- plot_func(safely_get_index(plot_obs_tm_nndm_ls, two_moons_lvl), 
                           plot_grid_tm_nndm, color_vals_tm,
                          size = 0.5, break_thresh = 0.08)
  
  
  pdc_nc_nndm <- plot_func(safely_get_index(plot_obs_nc_nndm_ls, circles_lvl), 
                           plot_grid_nc_nndm, color_vals_nc,
                           size = 0.5, break_thresh = 0.02)
  
  pdc_ts_nndm <- plot_func_tsne(safely_get_index(plot_obs_ts_nndm_ls, tsne_lvl), size=0.25)

  png(sprintf("output/toy_challenge/compare_ballet_clusterings-%s.png", file_suffix), 
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
    top = textGrob("BALLET Clustering Point Estimates", 
                   gp = gpar(fontsize = 18)),
    ncol = 3
  )
  dev.off()
}


plot_across_density_models(1,1,1, file_suffix="high")
plot_across_density_models(1,1,2, file_suffix="medium")
plot_across_density_models(2,2,3, file_suffix="low")
plot_across_density_models(3,3,4, file_suffix="elbow")

# Show the "elbow" figures for various data-set/models.

elbow_plot <- function(Ef, noise_frac) {
  d <- data.frame(density=Ef, log.density=log(Ef), ranks=rank(Ef))
  g <- ggplot(d, aes(x=ranks, y=log.density)) + geom_point(size=0.5) + 
    geom_vline(xintercept=noise_frac*length(Ef), color='red') +
    guides(color = 'none') + labs(x=NULL, y=NULL)
}

find_elbow_noise_frac <- function(plot_obs) {
  lvls <- str_match(colnames(plot_obs), 'db_pe_(.+)$')[,2]
  lvls <- lvls[!is.na(lvls)]
  
  # Analyze levels to decide on their ordering
  lvls_num <- as.numeric(lvls)
  # preset levels are in multiples of 0.05,
  # while an elbow level will not typically be..
  elbow_level <- which(floor(lvls_num*1000) %% 50 != 0)
  
  if(is.na(elbow_level)) {
    warning("Could not find elbow level in [", paste(lvls, collapse = ","),"]")
  } else {
    elbow_col <- paste0('db_pe_', lvls[elbow_level])
    attr(plot_obs[[elbow_col]], 'noise_frac')
  }
}


methods <- c("apt","dpmm", "nndm")
datasets <- c("two_moons", "circles", "tsne")

gr_obs <- rep(list(rep(list(NULL), length(datasets))), length(methods))

for (method in methods) {
  for (dataset in datasets) {
    
    sprintf("output/toy_challenge/density_pe_%s_%s.rds", method, dataset) |>
      readRDS() -> Ef
    
    sprintf("output/toy_challenge/plot_obs_%s_%s.rds", method, dataset) |>
      readRDS() |> find_elbow_noise_frac() -> noise_frac
    
    gr_obs[[method]][[dataset]] <- elbow_plot(Ef, noise_frac)
  }
}

gr_obs$dpmm$two_moons

#See https://cran.r-project.org/web/packages/gridExtra/vignettes/arrangeGrob.html

png("output/toy_challenge/ballet_elbow_plots.png", 
    width = 12, height = 8, units = 'in', res = 300)
grid.arrange(
  arrangeGrob(gr_obs$dpmm$two_moons,
              top = textGrob(
                "DP Mixture of Gaussians",
                gp = gpar(fontface = 3, fontsize = 14)
              ),
              left = textGrob("Two Moons", 
                          gp = gpar(fontface = 3, fontsize = 14), rot=90)),
  arrangeGrob(gr_obs$apt$two_moons,
              top = textGrob(
                "Adaptive Polya Tree",
                gp = gpar(fontface = 3, fontsize = 14)
              )),
  arrangeGrob(gr_obs$nndm$two_moons,
              top = textGrob(
                "NN Dirichlet Mixture",
                gp = gpar(fontface = 3, fontsize = 14)
              )),
  arrangeGrob(gr_obs$dpmm$circles,
              left=textGrob("Noisy Circles", gp = gpar(fontface = 3, fontsize = 14), 
                            rot=90)),
  gr_obs$apt$circles,
  gr_obs$nndm$circles,
  arrangeGrob(gr_obs$dpmm$tsne, left=textGrob("t-SNE", gp = gpar(fontface = 3, fontsize = 14), 
                                              rot=90)),
  gr_obs$apt$tsne, 
  gr_obs$nndm$tsne, 
  top = textGrob("Determining cutoff level using elbow plots", 
                 gp = gpar(fontsize = 18)),
  right = textGrob(TeX("log(\\hat{f}($x_i$))"), gp = gpar(fontsize = 14), rot=90),
  bottom = textGrob(TeX("rank of log(\\hat{f}($x_i$)) across $i \\in \\{1, \\ldots, n\\}$"), gp = gpar(fontsize = 14)),
  nrow = 3,
  ncol = 3
)
dev.off()


## Plot the various values for the tSNE dataset and the corresponding 
## persistent clustering..

plot_tnse_persistent <- function(plot_obs_ls) {
  plot_obs_ls |>
    discard_at(length(plot_obs_ls)) |> # remove the last element (the elbow)
      map("db_pe") |> rev() -> ctree_ls

  do.call(cbind, ctree_ls) -> ctree
  colnames(ctree) <- paste0("NoiseFrac", map_dbl(ctree_ls, \(x) attr(x,'noise_frac')))
  # In prep_labels, we replaced 0 to NAs. Change this back. 
  ctree[is.na(ctree)] <- 0
  clustree(ctree, prefix="NoiseFrac")

  select_persistent_clusters(ctree, 'NoiseFrac') |>
    rearrange_labels() -> pc
  pc <- factor(pc, levels=1:max(pc))
  attr(pc, 'noise_frac') <- NA #To denote persistent clustering

  #plot_obs_ls |>
  #  discard_at(length(plot_obs_ls)) |> 
  #    map(plot_func_tsne) -> gr_obs

  plot_obs_per <- plot_obs_ls[[1]]
  plot_obs_per$db_pe <- pc
  g_per <- plot_func_tsne(plot_obs_per)
  
  g_per
}

png(sprintf("output/toy_challenge/tsne-persistent-clustering.png"), 
    width = 12, height = 8, units = 'in', res = 300)
grid.arrange(
  arrangeGrob(plot_tnse_persistent(plot_obs_ts_ls),
    top = textGrob(
      "DP Mixture of Gaussians",
      gp = gpar(fontface = 3, fontsize = 14))
  ),
  arrangeGrob(plot_tnse_persistent(plot_obs_ts_apt_ls),
              top = textGrob(
                "Adaptive Polya Tree",
                gp = gpar(fontface = 3, fontsize = 14))
  ),
  arrangeGrob(plot_tnse_persistent(plot_obs_ts_nndm_ls),
              top = textGrob(
                "NN Dirichlet Mixture",
                gp = gpar(fontface = 3, fontsize = 14))
  ),
  ncol=2
)
dev.off()
