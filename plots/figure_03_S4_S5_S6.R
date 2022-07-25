library(dplyr)
library(ggplot2)
library(viridis)

# for plots, drop 2 elements of the list: Deepsea sole (doesn't converge) and
# Longspine thornyhead because it was run 2x

#---- plot results
null_index = readRDS("output/null_index_range15_sigma10.rds")
ll_index = readRDS("output/ll_index_range15_sigma10.rds")
#species_names = readRDS(file="species.rds")
#for(i in 1:length(null_index)) {
  #if(!is.null(null_index[[i]])) {
  #  null_index[[i]]$species <- species_names$common_name[i]
  #}
  #if(!is.null(ll_index[[i]])) {
  #  ll_index[[i]]$species <- species_names$common_name[i]
  #}
#}

null_df = bind_rows(null_index)
ll_df = bind_rows(ll_index)

null_df$model = "Constant"
ll_df$model = "Log-linear"
joined_df = rbind(ll_df, null_df)

joined_df = dplyr::filter(joined_df, common_name != "")

pdf("plots/Figure_S4_biomass_index_log.pdf")
ggplot(joined_df, aes(year, log_est, fill=model, col=model, group=model)) +
  geom_line(alpha=0.5) +
  geom_ribbon(aes(ymin = log(lwr), ymax = log(upr)), alpha=0.4, colour = NA) +
  facet_wrap(~ common_name, scale="free_y", ncol = 5) +
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
  xlab("") +
  scale_colour_brewer(palette="Dark2") +
  scale_fill_brewer(palette="Dark2")
  #scale_fill_viridis(discrete=TRUE,end=0.8) +
  #scale_color_viridis(discrete=TRUE,end=0.8)
dev.off()

joined_df$model = as.factor(joined_df$model)

top5_bottom5 = readRDS("output/top5_bottom5.rds")

# order based on species
sub_df <- joined_df
sub_df$common_name[which(sub_df$common_name == "Splitnose rockfish ")] = "Splitnose rockfish"
sub_df$common_name[which(sub_df$common_name == "Shortspine thornyhead ")] = "Shortspine thornyhead"

sub_df = dplyr::filter(sub_df, common_name %in% top5_bottom5$common_name)

sub_df$new_name = factor(sub_df$common_name, levels = c("Splitnose rockfish","Lingcod","Petrale sole","Darkblotched rockfish",
                                                     "Shortspine thornyhead"
                                                     ,"Rex sole","Spotted ratfish","Rosethorn rockfish","Longnose skate","Sandpaper skate"))

pdf("plots/Figure_3_biomass_index_top5.pdf",height = 5,width = 7)
ggplot(sub_df, aes(year, log_est, fill=model, col=model, group=model)) +
  geom_line(alpha=0.5) +
  geom_ribbon(aes(ymin = log(lwr), ymax = log(upr)), alpha=0.4, colour = NA) +
  facet_wrap(~ new_name, scale="free_y",nrow=2) +
  theme_bw() +
  ylab("Ln biomass index (+/- 2SE)") +
  theme(strip.background =element_rect(fill="white")) +
  theme(strip.text.x = element_text(size = 6),
        axis.text.x = element_text(size=5, angle = 90),
        axis.text.y = element_text(size=5),
        legend.justification = c(1, 0),
        legend.key.size = unit(0.3, "cm"),
        legend.text=element_text(size=rel(0.7)),
        legend.title=element_text(size=rel(0.8)),
        legend.position = c(0.99, 0.3)) +
  xlab("") +
  scale_colour_brewer(palette="Dark2") +
  scale_fill_brewer(palette="Dark2")
dev.off()

df = null_df[,c("year","est","common_name","se")]
df$est_ll = ll_df$est
df$se_ll = ll_df$se
df$ratio = df$est_ll / df$est
df$ratio_se = df$se_ll / df$se
df = dplyr::filter(df, common_name!="")

# join in names
# new_names = read.csv("data/name_change.csv") %>%
#   dplyr::rename(species = old_name)
# df = dplyr::left_join(df, new_names)
# df$new_name[which(is.na(df$new_name))] = "Shortspine thornyhead"

pdf("plots/Figure_S5_ratio.pdf")
ggplot(df, aes(year, ratio,group=common_name)) +
  geom_line() +
  facet_wrap(~ common_name, scale="free_y", ncol = 5) +
  theme_bw() +
  ylab("Ratio log-linear estimate / null estimated biomass") +
  theme_bw() +
  theme(strip.background =element_rect(fill="white")) +
  theme(strip.text.x = element_text(size = 6),
        axis.text.x = element_text(size=5, angle = 90),
        axis.text.y = element_text(size=5)) +
  xlab("")
dev.off()


pdf("plots/Figure_S6_ratio_se.pdf")
ggplot(df, aes(year, ratio_se,group=common_name)) +
  geom_line() +
  facet_wrap(~ common_name, scale="free_y", ncol = 5) +
  theme_bw() +
  ylab("Ratio of log-linear SE / null biomass SE") +
  theme_bw() +
  theme(strip.background =element_rect(fill="white")) +
  theme(strip.text.x = element_text(size = 6),
        axis.text.x = element_text(size=5, angle = 90),
        axis.text.y = element_text(size=5)) +
  xlab("")
dev.off()
