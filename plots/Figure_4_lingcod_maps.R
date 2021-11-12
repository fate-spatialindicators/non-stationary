library(sdmTMB)
library(dplyr)
library(ggplot2)
library(viridis)

# Species of interest
species = read.csv("survey_data/species_list.csv", fileEncoding="UTF-8-BOM")
names(species) = tolower(names(species))
species = dplyr::rename(species,
                        common_name = common.name,
                        scientific_name = scientific.name)

grid = readRDS("data/wc_grid.rds")
grid = dplyr::rename(grid, lon = X, lat = Y)
grid = dplyr::mutate(grid,
                     lon = lon*10, # scale to units of km
                     lat = lat*10,
                     depth_scaled = as.numeric(scale(-depth)),
                     depth_scaled2 = depth_scaled^2) %>%
  dplyr::select(-log_depth_scaled,
                -log_depth_scaled2)

grid$cell = seq(1,nrow(grid))
pred_grid = expand.grid(cell = grid$cell, year = 2003:2018)
pred_grid = dplyr::left_join(pred_grid, grid)
pred_grid$year = as.factor(pred_grid$year)
pred_grid$time = as.numeric(pred_grid$year) - floor(mean(unique(as.numeric(pred_grid$year))))

null_predictions = readRDS("output/null_predictions_summary.rds")
ll_predictions = readRDS("output/ll_predictions_summary.rds")


indx = which(species$species=="Lingcod")
# get predictions for each model
pred_grid$pred_ll_mean = ll_predictions[[indx]][,1]
pred_grid$pred_null_mean = null_predictions[[indx]][,1]
pred_grid$pred_ll_se = ll_predictions[[indx]][,2]
pred_grid$pred_null_se = null_predictions[[indx]][,2]
# subset years to 2003 and 2018
pred_grid = dplyr::filter(pred_grid, year %in%c(2003,2018))

pred_grid = dplyr::group_by(pred_grid, year) %>%
  dplyr::mutate(pred_ll_mean = (pred_ll_mean) - mean(pred_ll_mean),
                pred_null_mean = (pred_null_mean) - mean(pred_null_mean))

# 2 data frames, one for diff in means, one for diff in ses
diff_mean = dplyr::group_by(pred_grid, year, cell) %>%
  dplyr::summarise(diff = exp(pred_ll_mean - pred_null_mean),
                   lat = lat, lon = lon, metric="mean")
diff_se = dplyr::group_by(pred_grid, year, cell) %>%
  dplyr::summarise(diff = pred_ll_se / pred_null_se,
                   lat = lat, lon = lon, metric="se")

# make plots of difference
g1 = diff_mean %>%
  ggplot(aes(lon,lat,fill=diff)) +
  geom_tile() +
  scale_fill_gradient2(low="red",high="blue",midpoint=1, name="Ratio of means") +
  #scale_fill_continuous(name="Ratio of means") +
  #scale_fill_gradient2(high="blue",midpoint=0, name="Ratio of means") +
  facet_wrap(~year) +
  coord_fixed() +
  theme_bw() +
  xlab("") + ylab("") +
  theme(strip.background =element_rect(fill="white"))

g2 = diff_se %>%
  ggplot(aes(lon,lat,fill=(diff))) +
  geom_tile() +
  scale_fill_gradient2(low="red",high="blue",midpoint=1, name="Ratio of SEs") +
  facet_wrap(~year) +
  coord_fixed() +
  theme_bw() +
  xlab("") + ylab("") +
  theme(strip.background =element_rect(fill="white"))

library(grid)
library(gridExtra)
library(cowplot)
library(egg)
library(ggpubr)
#g0 = cowplot::plot_grid(g1,g2,nrow=2)
g0=egg::ggarrange(g1,g2,nrow=2)
pdf("plots/Figure_4_mean_cutoff_b.pdf", height = 7, width=7)
h0 = annotate_figure(g0,
                     left = textGrob("Latitude", rot = 90, vjust = 13, gp = gpar(cex = 1)),
                     bottom = textGrob("Longitude", hjust=1,vjust = -1,gp = gpar(cex =1)))
h0
dev.off()
