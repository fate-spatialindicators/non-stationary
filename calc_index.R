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
    null_index[[i]] <- get_index(null_predictions[[i]], bias_correct = TRUE)

    ll_predictions[[i]] <- predict(m_adult_ll, newdata = pred_grid, return_tmb_object = TRUE)
    ll_index[[i]] <- get_index(ll_predictions[[i]], bias_correct = TRUE)
  }
}
saveRDS(null_predictions,"output/null_predictions.rds")
saveRDS(ll_predictions,"output/ll_predictions.rds")
saveRDS(null_index,"output/null_index.rds")
saveRDS(ll_index,"output/ll_index.rds")

null_index = readRDS("output/null_index.rds")
ll_index = readRDS("output/ll_index.rds")

null_df = bind_rows(null_index[-8])
null_df$species = c(t(replicate(species$common_name[-c(8)],n=length(2003:2018))))

null_df$model = "Constant"
ll_df = bind_rows(ll_index[-8])
ll_df$species = c(t(replicate(species$common_name[-c(8)],n=length(2003:2018))))
ll_df$model = "Log-linear"
joined_df = rbind(ll_df, null_df)

pdf("plots/adult_biomass_index.pdf")
ggplot(joined_df, aes(year, log_est, color=model, group=model)) +
  geom_line() +
  geom_pointrange(aes(ymin=log_est-2*se, ymax = log_est+2*se), alpha=0.7,
    position = position_dodge(width = 0.5)) +
  facet_wrap(~ species, scale="free_y") +
  theme_bw() +
  scale_color_viridis(discrete=TRUE,end=0.8) +
  ylab("Log est (+/- 2SE)")
dev.off()
