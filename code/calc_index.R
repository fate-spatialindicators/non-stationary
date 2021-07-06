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

null_predictions = list()
ll_predictions = list()
null_index = list()
ll_index = list()

for(i in 1:nrow(species)){

  load(file=paste0("output/", sub(" ", "_", species$common_name[[i]]),"_ar1.RData"))

  est_index = TRUE
  if(class(ad_fit)=="try-error") est_index = FALSE
  if(class(ad_fit_ll)=="try-error") est_index = FALSE

  if(est_index==TRUE) {
    # calculate the biomass trend for the adult models
    null_predictions[[i]] <- predict(ad_fit, newdata = pred_grid, sims = 500)
    null_index[[i]] <- get_index_sims(null_predictions[[i]])

    ll_predictions[[i]] <- predict(ad_fit_ll, newdata = pred_grid, sims = 500)
    ll_index[[i]] <- get_index_sims(ll_predictions[[i]])
  }
  print(paste0("species ", i, " of ", nrow(species), " complete"))
}

#saveRDS(null_predictions,"output/null_predictions.rds") # too large to save
#saveRDS(ll_predictions,"output/ll_predictions.rds") # too large to save
saveRDS(null_index,"output/null_index.rds")
saveRDS(ll_index,"output/ll_index.rds")


#---- plot results
null_index = readRDS("output/null_index.rds")
ll_index = readRDS("output/ll_index.rds")

null_df = bind_rows(null_index)
null_df$species = c(t(replicate(species$common_name,n=length(2003:2018))))

null_df$model = "Constant"
ll_df = bind_rows(ll_index)
ll_df$species = c(t(replicate(species$common_name,n=length(2003:2018))))
ll_df$model = "Log-linear"
joined_df = rbind(ll_df, null_df)

joined_df = dplyr::filter(joined_df, species != "")
pdf("plots/biomass_index_log.pdf")
ggplot(joined_df, aes(year, log_est, fill=model, col=model, group=model)) +
  geom_line() +
  geom_ribbon(aes(ymin = log_est-se, ymax = log_est+se), alpha=0.4, colour = NA) +
  facet_wrap(~ species, scale="free_y", ncol = 5) +
  theme_bw() +
  ylab("Log estimate (+/- 1SE)") +
  theme(strip.background =element_rect(fill="white")) +
  theme(strip.text.x = element_text(size = 6),
        axis.text.x = element_text(size=5, angle = 90),
        axis.text.y = element_text(size=5),
        legend.position = c(1, 0),
        legend.justification = c(1, 0),
        legend.key.size = unit(0.3, "cm"),
        legend.text=element_text(size=rel(0.7)),
        legend.title=element_text(size=rel(0.8)))
dev.off()

joined_df$model = as.factor(joined_df$model)
pdf("plots/biomass_index_normal.pdf")
ggplot(joined_df, aes(year, exp(log_est), fill=model,group=model,col=model)) +
  geom_line() +
  geom_ribbon(aes(ymin = exp(log_est-se), ymax = exp(log_est+se)), alpha = 0.4, colour = NA) +
  facet_wrap(~ species, scale="free_y", ncol = 5) +
  theme_bw() +
  ylab("Estimate (+/- 1SE)") +
  theme(strip.background =element_rect(fill="white")) +
  theme(strip.text.x = element_text(size = 6),
        axis.text.x = element_text(size=5, angle = 90),
        axis.text.y = element_text(size=5),
        legend.position = c(1, 0),
        legend.justification = c(1, 0),
        legend.key.size = unit(0.3, "cm"),
        legend.text=element_text(size=rel(0.7)),
        legend.title=element_text(size=rel(0.8)))
dev.off()

df = null_df[,c("year","est","species")]
df$est_ll = ll_df$est
df$ratio = df$est_ll / df$est
df = dplyr::filter(df, species!="")

pdf("plots/ratio.pdf")
ggplot(df, aes(year, ratio,group=species)) +
  geom_line() +
  facet_wrap(~ species, scale="free_y", ncol = 5) +
  theme_bw() +
  ylab("Ratio log-linear estimate / null estimated biomass") +
  theme_bw() +
  theme(strip.background =element_rect(fill="white")) +
  theme(strip.text.x = element_text(size = 6),
        axis.text.x = element_text(size=5, angle = 90),
        axis.text.y = element_text(size=5))
dev.off()
