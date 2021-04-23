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
files <- paste0("output/", .sp, "_ar1.RData")
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
# .sp <- .sp[c(10, 14, 19)]
# files <- files[c(10, 14, 19)]
pred <- furrr::future_map2_dfr(.sp, files, .f = make_predictions)
# pred <- purrr::map2_dfr(.sp, files, .f = make_predictions)

pred %>%
  filter(year %in% c(2003L, 2018L)) %>%
  saveRDS(file = "output/predictions-2003-2018-ar1.rds")

plan(sequential) # avoid crashes

# x <- filter(pred, species == .sp[2], year %in% c(2003L, 2018L))
# ggplot(x, aes(lon, lat, fill = epsilon_st)) +
#   geom_raster() +
#   facet_grid(year~type) +
#   coord_fixed() +
#   scale_fill_gradient2()
# # ggsave("plots/dover-quadrant-maps-eps.png", width = 6, height = 8)
#
# ggplot(x, aes(lon, lat, fill = est)) +
#   geom_raster() +
#   facet_grid(year~type) +
#   coord_fixed() +
#   scale_fill_viridis_c()
# # ggsave("plots/dover-quadrant-maps-pred.png", width = 6, height = 8)
#
# x %>% group_by(lon, lat, year) %>%
#   summarise(diff = exp(est[type == "Null"]) / exp(est[type == "Non-stationary"])) %>%
#   ggplot(aes(lon, lat, fill = log10(diff))) +
#   geom_raster() +
#   facet_grid(~year) +
#   coord_fixed() +
#   scale_fill_gradient2()
#
# # g <- ggplot(x, aes(lon, lat, fill = omega_s)) +
# #   geom_raster() +
# #   facet_grid(year~type) +
# #   coord_fixed() +
# #   scale_fill_gradient2() +
# #   ggsidekick::theme_sleek()
# # ggsave("plots/dover-quadrant-maps-omega.png", width = 6, height = 8)
#
# x <- filter(pred, species == .sp[2], year %in% c(2003L, 2018L))
# x %>% group_by(year, type) %>%
#   filter(est > quantile(est, probs = 0.95)) %>%
#   ggplot(aes(lon, lat)) +
#   geom_raster(fill = "darkred") +
#   facet_grid(year~type) +
#   coord_fixed()

# x <- filter(pred, species == .sp[3], year %in% c(2003L, 2018L))
# x %>% group_by(year, type) %>%
#   filter(est > quantile(est, probs = 0.95)) %>%
#   ggplot(aes(lon, lat)) +
#   geom_raster(fill = "darkred") +
#   facet_grid(year~type) +
#   coord_fixed()

x <- filter(pred, year %in% c(2003L, 2018L))
dat <- x %>% group_by(year, type, species) %>%
  filter(est > quantile(est, probs = 0.95))

# models <- furrr::future_map(files, function(f) {
#   load(f)
#   ad_fit_ll
# })
#
# get_b_eps <- function(x) {
#   est <- as.list(x$sd_report, "Estimate", report = TRUE)$b_epsilon_logit
#   se <- as.list(x$sd_report, "Std. Error", report = TRUE)$b_epsilon_logit
#   tibble(
#     est_b_eps = est,
#     lwr_b_eps = est - 1.96 * se,
#     upr_b_eps = est + 1.96 * se
#   )
# }
#
# names(models) <- species$common_name
# b_eps <- furrr::future_map_dfr(models, get_b_eps, .id = "species")
# b_eps$species <- gsub(" ", "_", b_eps$species)
#
# dat <- left_join(dat, b_eps)

cols <- RColorBrewer::brewer.pal(5, "Dark2")[c(1, 4)]
# dir.create("plots/hotspots/", showWarnings = FALSE)
# for (i in unique(dat$species)) {
g <- dat %>%
  ggplot(aes(lon, lat, fill = type)) +
  geom_tile(alpha = 0.4) +
  facet_grid(species~year) +
  coord_fixed() +
  scale_fill_manual(values = cols) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank())
# ggtitle(i)
# g
ggsave(paste0("plots/hotspots-all.pdf"), width = 4, height = 80, limitsize = FALSE)
# }
