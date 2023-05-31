
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
  top = textGrob("Model-Based Clustering - W&G Credible Balls", 
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
               breaks = 0.07, col = 'black',
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
               breaks = 0.07, col = 'black',
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
