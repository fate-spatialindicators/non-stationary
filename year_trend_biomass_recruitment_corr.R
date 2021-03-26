
library(sdmTMB)
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)

# Species of interest
species_assmt = read.csv("data/Assessment_Time_Series_B_R.csv")
names(species_assmt) = tolower(names(species_assmt))
common_name = unique(species_assmt$species)
  
species_df = as.data.frame(common_name)
species_df$year_trend_adult = NA
species_df$year_se_adult = NA
species_df$year_trend_juv = NA
species_df$year_se_juv = NA
species_df$year_trend_total = NA
species_df$year_se_total = NA
species_df$year_trend_b = NA
species_df$year_trend_r = NA
species_df = filter(species_df, !common_name %in% c("Pacific hake")) # problem with old model...

for(j in 1:nrow(species_df)){
  comm_name = species_df$common_name[j]
  load(file=paste0("output/juv_ad_total_fits_all_models/", sub(" ", "_", comm_name),"_all_models.RData"))
  
  #indx = which(names(m_adult_ll$sd_report$value)=="b_epsilon")[1]
  #species_df$year_trend_adult[j] = m_adult_ll$sd_report$value[indx]
  #species_df$year_se_adultm_adult_ll] = m_adult_ll$sd_report$sd[indx]
  
  #indx = which(names(m_juv_ll$sd_report$value)=="b_epsilon")[1]
  #species_df$year_trend_juv[j] = m_juv_ll$sd_report$value[indx]
  #species_df$year_se_juv[j] = m_juv_ll$sd_report$sd[indx]
    
  indx = which(names(m_total_ll$sd_report$value)=="b_epsilon")[1]
  species_df[j, 'year_trend_total'] <- as.numeric(m_total_ll$sd_report$value[indx])
  species_df[j, 'year_se_total'] <- as.numeric(m_total_ll$sd_report$sd[indx])
  
  assmt = species_assmt %>% filter(species == comm_name) %>% mutate(b_b0 = b/max(b, na.rm=TRUE))
  species_df[j, 'year_trend_b'] = lm(b_b0 ~ year, data = assmt)$coefficients['year']
  #species_df[j, 'year_trend_r'] = lm(r_r0 ~ year, data = assmt)$coefficients['year']
}
species_df

p <- species_df %>%
  ggplot(aes(year_trend_b,year_trend_total)) +
  geom_point() +
  stat_smooth(method = lm) +
  xlab("Temporal trend in B/B0") +
  ylab("Temporal trend in spatiotemporal sd") +
  theme_bw() 
p  + ggpubr::stat_cor(method = "pearson", label.y = 0.1)

pdf("plots/figure_trend_epsilon_b.pdf")
p
dev.off()

jpeg("plots/figure_trend_epsilon_b.jpeg")
p
dev.off()