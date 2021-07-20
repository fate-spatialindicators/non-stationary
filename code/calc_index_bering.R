library(sdmTMB)
library(dplyr)
library(ggplot2)
library(viridis)
library(sp)
library(rgdal)

grid <- read.csv("data/ebs_grid.csv")
names(grid) <- tolower(names(grid))
grid <- dplyr::select(grid, lon, lat, shape_area)
coordinates(grid) <- c("lon", "lat")
proj4string(grid) <- CRS("+proj=longlat +datum=WGS84")
grid <- spTransform(grid, CRS("+proj=utm +zone=3"))
grid <- as.data.frame(grid)
grid$lon <- grid$lon/1000 # scale units to km
grid$lat <- grid$lat/1000

grid$cell <- seq(1, nrow(grid))
pred_grid <- expand.grid(cell = grid$cell, year = seq(1982, 2019))
pred_grid <- left_join(pred_grid, grid)
pred_grid$year <- as.factor(pred_grid$year)
pred_grid$time = as.numeric(pred_grid$year) - floor(mean(unique(as.numeric(pred_grid$year))))

load(file = "output/pollock_bering_ar1.RData")

null_predictions <- predict(ad_fit, newdata = pred_grid, sims = 500)
null_index <- get_index_sims(null_predictions)

ll_predictions <- predict(ad_fit_ll, newdata = pred_grid, sims = 500)
ll_index <- get_index_sims(ll_predictions)

saveRDS(null_predictions,"output/null_predictions_bering_pollock.rds")
saveRDS(ll_predictions,"output/ll_predictions_bering_pollock.rds")
saveRDS(null_index,"output/null_index_bering_pollock.rds")
saveRDS(ll_index,"output/ll_index_bering_pollock.rds")
rm(null_predictions, ll_predictions)

#---- plot results
#null_index = readRDS("output/null_index.rds")
#ll_index = readRDS("output/ll_index.rds")

null_index$model = "Constant"
ll_index$model = "Log-linear"
joined_df = rbind(ll_index, null_index)

pdf("plots/biomass_index_log_bering_pollock.pdf")
ggplot(joined_df, aes(year, log_est, fill=model, col=model, group=model)) +
  geom_line() +
  geom_ribbon(aes(ymin = log_est-se, ymax = log_est+se), alpha=0.4, colour = NA) +
  theme_bw() +
  ylab("Log estimate (+/- 1SE)") +
  theme(strip.background =element_rect(fill="white")) +
  theme(strip.text.x = element_text(size = 6),
        axis.text.x = element_text(size=5, angle = 90),
        axis.text.y = element_text(size=5))
dev.off()

joined_df$model = as.factor(joined_df$model)
pdf("plots/biomass_index_normal_bering_pollock.pdf")
ggplot(joined_df, aes(year, exp(log_est), fill=model,group=model,col=model)) +
  geom_line() +
  geom_ribbon(aes(ymin = exp(log_est-se), ymax = exp(log_est+se)), alpha = 0.4, colour = NA) +
  theme_bw() +
  ylab("Estimate (+/- 1SE)") +
  theme(strip.background =element_rect(fill="white")) +
  theme(strip.text.x = element_text(size = 6),
        axis.text.x = element_text(size=5, angle = 90),
        axis.text.y = element_text(size=5))
dev.off()
