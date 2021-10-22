library(dplyr)
library(ggplot2)
library(viridis)

# for plots, drop 2 elements of the list: Deepsea sole (doesn't converge) and
# Longspine thornyhead because it was run 2x
spp_to_drop = c(5, 31)
# Species of interest
species = read.csv("survey_data/species_list.csv", fileEncoding="UTF-8-BOM")
names(species) = tolower(names(species))
species = dplyr::rename(species,
                        common_name = common.name,
                        scientific_name = scientific.name)
species = species[-c(spp_to_drop),]


#---- plot results
null_index = readRDS("output/null_index_range15_sigma10.rds")
ll_index = readRDS("output/ll_index_range15_sigma10.rds")
null_index = null_index[-spp_to_drop]
ll_index = ll_index[-spp_to_drop]

null_df = bind_rows(null_index)
null_df$species = c(t(replicate(species$common_name,n=length(2003:2018))))

null_df$model = "Constant"
ll_df = bind_rows(ll_index)
ll_df$species = c(t(replicate(species$common_name,n=length(2003:2018))))
ll_df$model = "Log-linear"
joined_df = rbind(ll_df, null_df)

joined_df = dplyr::filter(joined_df, species != "")

# join in names
new_names = read.csv("data/name_change.csv") %>%
  dplyr::rename(species = old_name)
joined_df = dplyr::left_join(joined_df, new_names)
joined_df$new_name[which(is.na(joined_df$new_name))] = "Shortspine thornyhead"

pdf("plots/Figure_03_biomass_index_log.pdf")
ggplot(joined_df, aes(year, log_est, fill=model, col=model, group=model)) +
  geom_line() +
  geom_ribbon(aes(ymin = log(lwr), ymax = log(upr)), alpha=0.4, colour = NA) +
  facet_wrap(~ new_name, scale="free_y", ncol = 5) +
  theme_bw() +
  ylab("Ln biomass index (+/- 2SE)") +
  theme(strip.background =element_rect(fill="white")) +
  theme(strip.text.x = element_text(size = 6),
        axis.text.x = element_text(size=5, angle = 90),
        axis.text.y = element_text(size=5),
        legend.position = c(1, 0),
        legend.justification = c(1, 0),
        legend.key.size = unit(0.3, "cm"),
        legend.text=element_text(size=rel(0.7)),
        legend.title=element_text(size=rel(0.8))) +
  xlab("")
dev.off()

joined_df$model = as.factor(joined_df$model)
pdf("plots/Figure_S3_biomass_index_normal.pdf")
ggplot(joined_df, aes(year, exp(log_est), fill=model,group=model,col=model)) +
  geom_line() +
  geom_ribbon(aes(ymin = exp(log_est-se), ymax = exp(log_est+se)), alpha = 0.4, colour = NA) +
  facet_wrap(~ new_name, scale="free_y", ncol = 5) +
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
        legend.title=element_text(size=rel(0.8))) +
  xlab("")
dev.off()

df = null_df[,c("year","est","species")]
df$est_ll = ll_df$est
df$ratio = df$est_ll / df$est
df = dplyr::filter(df, species!="")

# join in names
new_names = read.csv("data/name_change.csv") %>%
  dplyr::rename(species = old_name)
df = dplyr::left_join(df, new_names)
df$new_name[which(is.na(df$new_name))] = "Shortspine thornyhead"

pdf("plots/Figure_S4_ratio.pdf")
ggplot(df, aes(year, ratio,group=new_name)) +
  geom_line() +
  facet_wrap(~ new_name, scale="free_y", ncol = 5) +
  theme_bw() +
  ylab("Ratio log-linear estimate / null estimated biomass") +
  theme_bw() +
  theme(strip.background =element_rect(fill="white")) +
  theme(strip.text.x = element_text(size = 6),
        axis.text.x = element_text(size=5, angle = 90),
        axis.text.y = element_text(size=5)) +
  xlab("")
dev.off()
