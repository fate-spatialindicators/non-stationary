# WORK IN PROGRESS! 2021-04-09 SA

library(sdmTMB)
library(dplyr)
library(ggplot2)

species <- read.csv("survey_data/species_list.csv")
names(species) <- tolower(names(species))
species <- rename(species,
  common_name = common.name,
  scientific_name = scientific.name
)

grid <- readRDS("data/wc_grid.rds")
grid <- rename(grid, lon = X, lat = Y)
grid <- mutate(grid,
  depth_scaled = as.numeric(scale(-depth)),
  depth_scaled2 = depth_scaled^2
) %>%
  select(
    -log_depth_scaled,
    -log_depth_scaled2
  )
grid$cell <- seq(1, nrow(grid))
pred_grid <- expand.grid(cell = grid$cell, year = seq(2003L, 2018L))
pred_grid <- left_join(pred_grid, grid, by = "cell")
pred_grid$year <- as.factor(pred_grid$year)

null_predictions <- list()
ll_predictions <- list()

for (i in seq(1, nrow(species))) {
  comm_name <- species$common_name[i]
  load(file = paste0("output/",
    sub(" ", "_", comm_name), "_all_models.RData"))
  null_predictions[[i]] <- predict(ad_fit, newdata = pred_grid)
  ll_predictions[[i]] <- predict(ad_fit_ll, newdata = pred_grid)
}
saveRDS(null_predictions,"output/null_predictions_hotspot.rds")
saveRDS(ll_predictions,"output/ll_predictions_hotsplot.rds")

null_predictions <- readRDS("output/null_predictions.rds")
ll_predictions <- readRDS("output/ll_predictions.rds")

# look at predictions for splitnose rockfish, which are a case where there's a difference in biomass index
null_pred <- null_predictions[[13]]$data
ll_pred <- ll_predictions[[13]]$data
ll_pred$null_est <- null_pred$est
ll_pred$ll_est <- ll_pred$est

p1 <- ggplot(ll_pred, aes(null_est, ll_est)) +
  geom_point(col = "darkblue", alpha = 0.3) +
  geom_abline(aes(intercept = 0, slope = 1), col = "grey", alpha = 0.6) +
  facet_wrap(~year) +
  xlab("Log density, null model") +
  ylab("Log density, Log-linear model") +
  theme_bw()

pdf("plots/null_vs_ll_estimates_splitnose_rockfish.pdf")
p1
dev.off()
