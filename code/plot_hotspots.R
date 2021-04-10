# WORK IN PROGRESS! 2021-04-09 SA

library(sdmTMB)
library(dplyr)
library(ggplot2)
library(future)
plan(multisession, workers = max(floor(availableCores() / 4), 4L))

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
)
grid$cell <- seq(1, nrow(grid))
pred_grid <- expand.grid(cell = grid$cell, year = c(2003L, 2018L))
pred_grid <- left_join(pred_grid, grid, by = "cell")
pred_grid$year <- factor(pred_grid$year, levels = as.character(seq(2003, 2018)))

.sp <- sub(" ", "_", species$common_name)
files <- paste0("output/", .sp, "_all_models.RData")
make_predictions <- function(sp, f) {
  load(f)
  p1 <- predict(ad_fit, newdata = pred_grid)
  p2 <- predict(ad_fit_ll, newdata = pred_grid)
  p1$type <- "Null"
  p2$type <- "Non-stationary"
  bind_rows(p1, p2) %>%
    mutate(species = sp)
}
pred <- furrr::future_map2(.sp, files, .f = make_predictions)
