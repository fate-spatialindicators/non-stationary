# Species of interest
species_names = readRDS(file="species.rds")
grid = readRDS("data/wc_grid.rds")
grid = dplyr::rename(grid, lon = X, lat = Y)
grid = dplyr::mutate(grid,
                     lon = lon*10, # scale to units of km
                     lat = lat*10,
                     depth_scaled = as.numeric(scale(-depth)),
                     depth_scaled2 = depth_scaled^2) %>%
  dplyr::select(-log_depth_scaled,
                -log_depth_scaled2)

grid$cell = seq(1,nrow(grid))
pred_grid = expand.grid(cell = grid$cell, year = 2003:2018)
pred_grid = dplyr::left_join(pred_grid, grid)
pred_grid$fyear = as.factor(pred_grid$year)
pred_grid$time = as.numeric(as.character(pred_grid$year)) - floor(mean(unique(as.numeric(as.character(pred_grid$year)))))
#as.numeric(pred_grid$year) - mean(pred_grid$year)
pred_grid$year = as.factor(pred_grid$year)

null_predictions_summ = list()
ll_predictions_summ = list()
null_index = list()
ll_index = list()

for(ii in 28:nrow(species)){
print(ii)
  load(file=paste0("output/", sub(" ", "_", species$common_name[[ii]]), "_ar1_priors.RData"))

  df <- data.frame(resid = residuals(ad_fit_ll, type="mvn-laplace"),
                   common_name = ad_fit_ll$data$common_name[ii])
  if(ii==1) {
    all_df = df
  } else {
    all_df = rbind(all_df, df)
  }
}

pdf("plots/Figure_S8.pdf")
ggplot(all_df, aes(sample=resid)) +
  geom_qq(alpha=0.1) + geom_qq_line(col='blue') +
  facet_wrap(~common_name) +
  xlab("Theoretical") + ylab("Sample") +
  theme(strip.text.x = element_text(size =6))
dev.off()
