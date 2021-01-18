library(sdmTMB)
library(dplyr)
library(ggplot2)

# Species of interest
species = read.csv("survey_data/species_list_revised.csv")
names(species) = tolower(names(species))
species = dplyr::rename(species,
  common_name = common.name,
  scientific_name = scientific.name,
  juv_threshold = max.length.cm)

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

null_predictions = list()
ll_predictions = list()
null_index = list()
ll_index = list()

for(i in 1:nrow(species)){

  comm_name = species$common_name[i]
  dat = readRDS(paste0("data/", sub(" ", "_", comm_name), "_expanded.rds"))

  dat$depth_scaled = as.numeric(scale(dat$depth_m))
  dat$depth_scaled2 = dat$depth_scaled^2

  load(file=paste0("output/", sub(" ", "_", comm_name),"_all_models.RData"))

  est_index = TRUE
  if(class(m_adult)=="try-error") est_index = FALSE
  if(class(m_adult_ll)=="try-error") est_index = FALSE

  if(est_index==TRUE) {
    # calculate the biomass trend for the adult models
    null_predictions[[i]] <- predict(m_adult, newdata = pred_grid, return_tmb_object = TRUE)
    #null_index[[i]] <- get_index(null_predictions[[i]], bias_correct = TRUE)

    ll_predictions[[i]] <- predict(m_adult_ll, newdata = pred_grid, return_tmb_object = TRUE)
    #ll_index[[i]] <- get_index(ll_predictions[[i]], bias_correct = TRUE)
  }
}
#saveRDS(null_predictions,"output/null_predictions.rds")
#saveRDS(ll_predictions,"output/ll_predictions.rds")
#saveRDS(null_index,"output/null_index.rds")
#saveRDS(ll_index,"output/ll_index.rds")

null_predictions = readRDS("output/null_predictions.rds")
ll_predictions = readRDS("output/ll_predictions.rds")

# look at predictions for splitnose rockfish, which are a case where there's a difference in biomass index
null_pred = null_predictions[[13]]$data
ll_pred = ll_predictions[[13]]$data
ll_pred$null_est = null_pred$est
ll_pred$ll_est = ll_pred$est

p1 = ggplot(ll_pred, aes(null_est, ll_est)) +
  geom_point(col="darkblue", alpha = 0.3) +
  geom_abline(aes(intercept=0,slope=1),col="grey",alpha=0.6) +
  facet_wrap(~year) +
  xlab("Log density, null model") +
  ylab("Log density, Log-linear model") +
  theme_bw()

pdf("plots/null_vs_ll_estimates_splitnose_rockfish.pdf")
p1
dev.off()
