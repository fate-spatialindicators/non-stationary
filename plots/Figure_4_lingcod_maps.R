library(sdmTMB)
library(dplyr)
library(ggplot2)
library(viridis)
library(grid)
library(gridExtra)
library(cowplot)
library(egg)
library(ggpubr)

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

# commented out -- but the de-meaning is providing a geometric mean
# pred_grid = dplyr::group_by(pred_grid, year) %>%
#   dplyr::mutate(pred_ll_mean = (pred_ll_mean) - mean(pred_ll_mean),
#                 pred_null_mean = (pred_null_mean) - mean(pred_null_mean))
#
#

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
  facet_grid(~year,space="free") +
  coord_fixed() +
  theme_bw() +
  xlab("") + ylab("") +
  theme(strip.background =element_rect(fill="white"))

g2 = diff_se %>%
  ggplot(aes(lon,lat,fill=(diff))) +
  geom_tile() +
  scale_fill_gradient2(low="red",high="blue",midpoint=1, name="Ratio of SEs") +
  facet_grid(~year,space="free") +
  coord_fixed() +
  theme_bw() +
  xlab("") + ylab("") +
  theme(strip.background =element_rect(fill="white"))

#4448.492 is the cutoff at 40 deg 10 min (40.16667)
g3 = g1 + ylim(4448.492,5400) + xlim(300,450)
g4 = g2 + ylim(4448.492,5400) + xlim(300,450)


pdf("plots/Figure_S7.pdf")
gA <- ggplotGrob(g1)
gB <- ggplotGrob(g2)
h0 = annotate_figure(rbind(gA, gB),
                     left = textGrob("Northings", rot = 90, vjust = 13, gp = gpar(cex = 1)),
                     bottom = textGrob("Eastings", hjust=1,vjust = -1,gp = gpar(cex =1)))
h0
dev.off()

#g0 = cowplot::plot_grid(g1,g2,nrow=2)
#g0=egg::ggarrange(g1,g2,nrow=2)
g3 <- g3 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
g4 <- g4 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("plots/Figure_4.pdf")
gA <- ggplotGrob(g3)
gB <- ggplotGrob(g4)
h0 = annotate_figure(rbind(gA, gB),
                     left = textGrob("Northings", rot = 90, vjust = 13, gp = gpar(cex = 1)),
                     bottom = textGrob("Eastings", hjust=1,vjust = -1,gp = gpar(cex =1)))
h0
dev.off()
