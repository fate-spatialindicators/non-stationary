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
pred_grid$fyear = as.factor(pred_grid$year)
pred_grid$time = as.numeric(as.character(pred_grid$year)) - floor(mean(unique(as.numeric(as.character(pred_grid$year)))))
#as.numeric(pred_grid$year) - mean(pred_grid$year)
pred_grid$year = as.factor(pred_grid$year)

null_predictions_summ = list()
ll_predictions_summ = list()
null_index = list()
ll_index = list()

for(i in 1:nrow(species)){

  load(file=paste0("output/", sub(" ", "_", species$species[[i]]),"_ar1_priors.RData"))

  est_index = TRUE
  if(class(ad_fit)=="logical") est_index = FALSE
  if(class(ad_fit_ll)=="logical") est_index = FALSE

  if(est_index==TRUE) {
    null_predictions <- predict(ad_fit, newdata = pred_grid, nsim = 500)
    mean_null <- apply(null_predictions, 1, mean)
    sd_null <- apply(null_predictions, 1, sd)
    null_predictions_summ[[i]] <- cbind(mean_null, sd_null)
    null_index[[i]] <- get_index_sims(null_predictions)
    null_index[[i]]$common_name <- ad_fit$data$common_name[1]

    ll_predictions <- predict(ad_fit_ll, newdata = pred_grid, nsim = 500)
    mean_ll <- apply(ll_predictions, 1, mean)
    sd_ll <- apply(ll_predictions, 1, sd)
    ll_predictions_summ[[i]] <- cbind(mean_ll, sd_ll)
    ll_index[[i]] <- get_index_sims(ll_predictions)
    ll_index[[i]]$common_name <- ad_fit_ll$data$common_name[1]

    # also calculate epsilon_st
    pred_null <- predict(ad_fit, newdata = pred_grid)
    pred_null$common_name = species$common_name[[i]]
    pred_null$species = species$scientific_name[[i]]
    pred_null$model = "Null"
    pred_ll <- predict(ad_fit_ll, newdata = pred_grid)
    pred_ll$common_name = species$common_name[[i]]
    pred_ll$species = species$scientific_name[[i]]
    pred_ll$model = "Log-linear"
    if(i==1) {
      df_out = rbind(pred_null, pred_ll)
    } else {
      df_out = rbind(df_out, rbind(pred_null, pred_ll))
    }
  }

  print(paste0("species ", i, " of ", nrow(species), " complete"))
}

saveRDS(df_out, "output/predictions_all.rds")
saveRDS(null_predictions_summ,"output/null_predictions_summary.rds")
saveRDS(ll_predictions_summ,"output/ll_predictions_summary.rds")
saveRDS(null_index,"output/null_index_range15_sigma10.rds")
saveRDS(ll_index,"output/ll_index_range15_sigma10.rds")


