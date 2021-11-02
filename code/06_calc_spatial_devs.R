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

indx = which(species$species=="Darkblotched rockfish")

#for(i in 1:nrow(species)){
  # join predictions to data frame
  pred_grid$pred_ll = ll_predictions[[indx]][,1]
  pred_grid$pred_null = null_predictions[[indx]][,1]

  # calculate de-meaned data frame, subtracting off spatial mean for each cell
  pred_grid <- dplyr::group_by(pred_grid, cell) %>%
    dplyr::mutate(ll_mean = mean(pred_ll),
                  null_mean = mean(pred_null)) %>%
    dplyr::ungroup()
  pred_grid$ll_dev <- pred_grid$pred_ll - pred_grid$ll_mean
  pred_grid$null_dev <- pred_grid$pred_null - pred_grid$null_mean
  pred_grid = dplyr::select(pred_grid, -ll_mean, -null_mean)

  # subset years to 2003 and 2018
  pred_grid = dplyr::filter(pred_grid, year %in%c(2003,2018))

  pred_grid2 <- pred_grid
  pred_grid2$est = pred_grid2$null_dev
  pred_grid2$type = "Null"
  pred_grid$est = pred_grid$ll_dev
  pred_grid$type = "Non-stationary"
  pred_grid = rbind(pred_grid, pred_grid2)
  pred_grid$model = paste0(pred_grid$year, pred_grid$type)

  library(ggplot2)
# first version: this plots fields that are spatially de-meaned
  g0 =  ggplot(pred_grid, aes(lon,lat,fill=est)) +
    geom_tile() +
    scale_fill_gradient2() +
    facet_grid(year~type) +
    coord_fixed() +
    theme_bw()

  # second version: this subtracts off the mean for each year / model
  g1 =  dplyr::group_by(pred_grid, model, year) %>%
    dplyr::mutate(est = est - mean(est)) %>%
  ggplot(aes(lon,lat,fill=est)) +
    geom_tile() +
    scale_fill_gradient2() +
    facet_grid(year~type) +
    coord_fixed() +
    theme_bw()

  pdf("plots/Figure_4.pdf")
  g1
  dev.off()
