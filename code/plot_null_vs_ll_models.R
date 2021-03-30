
# Species of interest
species = read.csv("survey_data/species_list.csv")
names(species) = tolower(names(species))
species = dplyr::rename(species,
                        common_name = common.name,
                        scientific_name = scientific.name)

comm_name = "lingcod" # "rockfish, splitnose"

load(paste0("output/",comm_name,"_all_models.Rdata"))

# load models
fit = ad_fit
fit_ll = ad_fit_ll

# load full dataset
all_dat = readRDS("survey_data/all_data.rds")
sub = dplyr::filter(all_dat, scientific_name == species$scientific_name[i])
sub$depth_scaled = as.numeric(scale(sub$depth_m))
sub$depth_scaled2 = sub$depth_scaled^2

# create prediction grid
grid = readRDS("data/wc_grid.rds")
grid = dplyr::rename(grid, lon = X, lat = Y)
grid = dplyr::mutate(grid,
                     depth_scaled = as.numeric(scale(-depth)),
                     depth_scaled2 = depth_scaled^2) %>%
  dplyr::select(-log_depth_scaled,
                -log_depth_scaled2)

grid$cell = seq(1,nrow(grid))
pred_grid = expand.grid(cell = grid$cell, year = 2003:2018)
pred_grid = dplyr::left_join(pred_grid, grid)
pred_grid$lon = pred_grid$lon*10
pred_grid$lat = pred_grid$lat*10
pred_grid$year = as.factor(pred_grid$year)

# generate predictions
pred <- predict(ad_fit, newdata = pred_grid)
pred_ll <- predict(ad_fit_ll, newdata = pred_grid)
pred_diff = pred
pred_diff[,"est"] = pred[,"est"] - pred_ll[,"est"]
pred_diff[,"est_rf"] = pred[,"est_rf"] - pred_ll[,"est_rf"]
pred_diff[,"omega_s"] = pred[,"omega_s"] - pred_ll[,"omega_s"]
pred_diff[,"zeta_s"] = pred[,"zeta_s"] - pred[,"zeta_s"]
pred_diff[,"epsilon_st"] = pred[,"epsilon_st"] - pred[,"epsilon_st"]

# plot spatial component
p1 = dplyr::filter(pred, year==levels(pred$year)[1]) %>%
  ggplot(aes(lon, lat, fill = omega_s)) + geom_tile() +
  scale_fill_gradient2() + ggtitle("Spatial field, null model")
p2 = dplyr::filter(pred_diff, year==levels(pred$year)[1]) %>%
  ggplot(aes(lon, lat, fill = omega_s)) + geom_tile() +
  scale_fill_gradient2() + ggtitle("Diff, null spatial field - ll spatial field")

max_year = length(levels(pred$year)) - 1
# plot spatiotemporal field in first and last time step
pred_1 = dplyr::filter(pred, year%in%c(levels(pred$year)[1], levels(pred$year)[max_year]))
pred_2 = dplyr::filter(pred_ll, year%in%c(levels(pred$year)[1], levels(pred$year)[max_year]))
pred_1$model = "null"
pred_2$model = "loglinear"
pred_all = rbind(pred_1, pred_2)
pred_all$combo = paste0(pred_all$model, ": ", pred_all$year)

p3 = pred_all %>%
  ggplot(aes(lon, lat, fill = epsilon_st)) +
  geom_tile() +
  scale_fill_gradient2() +
  facet_wrap(~combo) +
  ggtitle("Estimated spatiotemporal field")

p4 = pred_all %>%
  ggplot(aes(lon, lat, fill = est_rf)) +
  geom_tile() +
  scale_fill_gradient2() +
  facet_wrap(~combo) +
  ggtitle("Estimated spatial + st field")

p5 = pred_all %>%
  ggplot(aes(lon, lat, fill = est_non_rf)) +
  geom_tile() +
  scale_fill_gradient2() +
  facet_wrap(~combo) + ggtitle("Estimate - non rf")

p6 = pred_all %>%
  ggplot(aes(lon, lat, fill = est)) +
  geom_tile() +
  scale_fill_gradient2() +
  facet_wrap(~combo) + ggtitle("Prediction (log space)")



pdf(paste0("plots/difference_null_ll_models_",comm_name,".pdf"))
gridExtra::grid.arrange(p1,p2,nrow=1)

p3
p4
p5
p6
dev.off()

