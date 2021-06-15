library(sdmTMB)
library(dplyr)
library(ggplot2)
library(viridis)

# Species of interest
species = read.csv("survey_data/species_list_subregional_indices.csv", fileEncoding="UTF-8-BOM")
names(species) = tolower(names(species))
species = dplyr::rename(species,
                        common_name = common.name,
                        scientific_name = scientific.name)

grid = readRDS("data/wc_grid.rds")
grid = dplyr::rename(grid, lon = X, lat = Y)
grid = dplyr::mutate(grid,
  lon = lon*10,
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

# create subsets of grid to predict to biogeographic regions
pred_grid_N = filter(pred_grid, lat > 4500)
pred_grid_S = filter(pred_grid, lat <= 4500)
#pred_grid_C = filter(pred_grid, lat <= 4500, lat >= 3850)


null_predictions_N = list()
ll_predictions_N = list()
null_index_N = list()
ll_index_N = list()
null_predictions_S = list()
ll_predictions_S = list()
null_index_S = list()
ll_index_S = list()

for(i in 1:nrow(species)){
  load(file=paste0("output/", sub(" ", "_", species$common_name[i]),"_ar1.RData"))

  est_index = TRUE
  if(class(ad_fit)=="try-error") est_index = FALSE
  if(class(ad_fit_ll)=="try-error") est_index = FALSE

  if(est_index==TRUE) {
    # calculate the biomass trend for the adult models, by subregion
    null_predictions_N[[i]] <- predict(ad_fit, newdata = pred_grid_N, return_tmb_object = TRUE)
    null_index_N[[i]] <- tryCatch({get_index(null_predictions_N[[i]], bias_correct = TRUE)}, error = function(e) NA)
    null_predictions_S[[i]] <- predict(ad_fit, newdata = pred_grid_S, return_tmb_object = TRUE)
    null_index_S[[i]] <- tryCatch({get_index(null_predictions_S[[i]], bias_correct = TRUE)}, error = function(e) NA)

    ll_predictions_N[[i]] <- predict(ad_fit_ll, newdata = pred_grid_N, return_tmb_object = TRUE)
    ll_index_N[[i]] <- tryCatch({get_index(ll_predictions_N[[i]], bias_correct = TRUE)}, error = function(e) NA)
    ll_predictions_S[[i]] <- predict(ad_fit_ll, newdata = pred_grid_S, return_tmb_object = TRUE)
    ll_index_S[[i]] <- tryCatch({get_index(ll_predictions_S[[i]], bias_correct = TRUE)}, error = function(e) NA)
  }
  print(paste0("species ", i, " of ", nrow(species), " complete"))
}
saveRDS(null_index_N,"output/null_index_N.rds")
saveRDS(ll_index_N,"output/ll_index_N.rds")
saveRDS(null_index_S,"output/null_index_S.rds")
saveRDS(ll_index_S,"output/ll_index_S.rds")


# null_index_N = readRDS("output/null_index_N.rds")
# ll_index_N = readRDS("output/ll_index_N.rds")
# null_index_S = readRDS("output/null_index_S.rds")
# ll_index_S = readRDS("output/ll_index_S.rds")

null_df = bind_rows(null_index_N, null_index_S)
null_df$species = rep(c(t(replicate(species$common_name, n=length(2003:2018)))), 2)
null_df$region = rep(c(t(replicate("N", n=length(2003:2018))), t(replicate("S", n=length(2003:2018)))), nrow(species))
null_df$model = "Constant"

ll_df = bind_rows(ll_index_N, ll_index_S)
ll_df$species = null_df$species
ll_df$region = null_df$region
ll_df$model = "Log-linear"

joined_df = rbind(ll_df, null_df)
joined_df = dplyr::filter(joined_df, species != "")
pdf("plots/biomass_index_log_subregional.pdf")
ggplot(joined_df, aes(year, log_est, fill=model, col=model, shape=region, group=interaction(model,region))) +
  geom_line() +
  geom_ribbon(aes(ymin=log_est-se, ymax = log_est+se), alpha=0.4) +
  facet_wrap(~ species, scale="free_y") +
  theme_bw() +
  ylab("Log estimate (+/- 1SE)") +
  theme(strip.background =element_rect(fill="white")) +
  theme(strip.text.x = element_text(size=4),
        axis.text.x = element_text(size=5, angle = 90))
dev.off()

joined_df$model = as.factor(joined_df$model)
pdf("plots/biomass_index_normal_subregional.pdf")
ggplot(joined_df, aes(year, exp(log_est), fill=model, col=model, shape=region, group=interaction(model,region))) +
  geom_line() +
  geom_ribbon(aes(ymin=exp(log_est-se), ymax = exp(log_est+se)), alpha=0.4) +
  facet_wrap(~ species, scale="free_y") +
  theme_bw() +
  ylab("Estimate (+/- 1SE)") +
  theme(strip.background =element_rect(fill="white")) +
  theme(strip.text.x = element_text(size=4),
        axis.text.x = element_text(size=5, angle = 90))
dev.off()

df = null_df[,c("year","est","species","region")]
df$est_ll = ll_df$est
df$ratio = df$est_ll / df$est
df = dplyr::filter(df, species!="")

pdf("plots/ratio_subregional.pdf")
ggplot(df, aes(year, ratio, col=region, group=interaction(model,region))) +
  geom_line() +
  facet_wrap(~ species, scale="free_y") +
  theme_bw() +
  ylab("Ratio log-linear estimate / null estimated biomass") +
  theme_bw() +
  theme(strip.background =element_rect(fill="white")) +
  theme(strip.text.x = element_text(size = 4),
        axis.text.x = element_text(size=5, angle = 90))
dev.off()


