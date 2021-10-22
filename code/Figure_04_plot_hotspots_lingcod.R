library(sdmTMB)
library(dplyr)
library(ggplot2)
library(future)
theme_set(ggsidekick::theme_sleek())
plan(multisession, workers = max(floor(availableCores() / 4), 4L))

species <- read.csv("survey_data/species_list.csv", fileEncoding="UTF-8-BOM")
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
pred_grid <- expand.grid(cell = grid$cell, year = seq(2003L, 2018L))
pred_grid <- left_join(pred_grid, grid, by = "cell")
pred_grid$year <- as.factor(pred_grid$year)

pred_grid$time <- as.numeric(pred_grid$year) -
  floor(mean(unique(as.numeric(pred_grid$year))))

pred_grid$lat <- pred_grid$lat * 10
pred_grid$lon <- pred_grid$lon * 10

.sp <- sub(" ", "_", species$common_name)
# files <- paste0("output/", .sp, "_all_models.RData")
files <- paste0("output/", .sp, "_ar1_priors.RData")
make_predictions <- function(sp, f) {
  cat(sp, "\n")
  load(f)
  # in case missing in old model:
  ad_fit$epsilon_predictor <- "time"
  ad_fit_ll$epsilon_predictor <- "time"

  p1 <- predict(ad_fit, newdata = pred_grid)
  p2 <- predict(ad_fit_ll, newdata = pred_grid)
  p1$type <- "Null"
  p2$type <- "Non-stationary"
  bind_rows(p1, p2) %>%
    mutate(species = sp)
}

pred <- make_predictions(.sp[10], files[10])
#pred <- furrr::future_map2_dfr(.sp, files, .f = make_predictions)
# pred <- purrr::map2_dfr(.sp, files, .f = make_predictions)

pred %>%
  filter(year %in% c(2003L, 2018L)) %>%
  saveRDS(file = "output/predictions-2003-2018-ar1.rds")

plan(sequential) # avoid crashes


x <- filter(pred, year %in% c(2003L, 2018L))
dat <- x %>% group_by(year, type, species) %>%
  filter(est > quantile(est, probs = 0.95))

cols <- RColorBrewer::brewer.pal(5, "Dark2")[c(1, 4)]
# dir.create("plots/hotspots/", showWarnings = FALSE)
# for (i in unique(dat$species)) {
g <- dat %>%
  ggplot(aes(lon, lat, fill = type)) +
  geom_tile(alpha = 0.1) +
  facet_grid(~year) +
  coord_fixed() +
  scale_fill_manual(values = cols) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank())
# ggtitle(i)
pdf("plots/Figure_04_hotspots-lingcod.pdf")
g
dev.off()
#ggsave(paste0(""), width = 4, height = 80, limitsize = FALSE)
# }
