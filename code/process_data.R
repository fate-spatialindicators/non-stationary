library(sdmTMB)
library(dplyr)
library(ggplot2)
library(viridis)

# Species of interest
species = read.csv("survey_data/species_list.csv")
names(species) = tolower(names(species))
species = dplyr::rename(species,
  common_name = common.name,
  scientific_name = scientific.name)

for(i in 1:nrow(species)){

  comm_name = species$common_name[i]

  load(file=paste0("output/", sub(" ", "_", comm_name),"_all_models.RData"))

  df = data.frame(name = comm_name,
    model = c("adult", "adult"),
    loglinear = c(FALSE,TRUE))
  df$aic = NA
  df$trend = NA
  df$trend_se = NA

  if(class(ad_fit)!="try-error") {
    df$aic[1] = AIC(ad_fit)
  }
  if(class(ad_fit_ll)!="try-error") {
    df$aic[2] = AIC(ad_fit_ll)
    df$trend[2] = ad_fit_ll$sd_report$value[which(names(ad_fit_ll$sd_report$value) == "b_epsilon_logit")]
    df$trend_se[2] = ad_fit_ll$sd_report$sd[which(names(ad_fit_ll$sd_report$value) == "b_epsilon_logit")]
  }

  if(i==1) {
    df_all = df
  } else {
    df_all = rbind(df_all, df)
  }
}


# look at coefficients for adult models
jpeg("plots/Trend_summaries.jpeg")

p1 = dplyr::filter(df_all, loglinear==TRUE, model=="adult", name %in% c("","thornyhead, longspine","sole, deepsea","skate, sandpaper")==FALSE) %>%
  ggplot(aes(name,trend)) +
  geom_pointrange(aes(ymin=trend-2*trend_se, ymax=trend+2*trend_se),col="darkblue",size=0.8) +
  xlab("Species") + ylab("Trend in adult spatiotemporal sd (+/- 2SE)") +
  geom_hline(aes(yintercept=0),col="red",alpha=0.6) +
  coord_flip() +
  theme_bw()
p1
dev.off()

df = dplyr::filter(df_all, !is.na(trend)) %>%
  dplyr::rename(common_name = name)
df = dplyr::left_join(species, df)

jpeg("plots/Trend_depletion.jpeg")

ggplot(dplyr::filter(df, common_name %in% c("","thornyhead, longspine","sole, deepsea")==FALSE), aes(depletion,trend)) +
  geom_pointrange(aes(ymin=trend-trend_se,ymax=trend+trend_se)) +
  #geom_point() +
  geom_smooth(method="lm") +
  xlim(0.25,1) +
  ylim(-0.25,0.11)
dev.off()

p2 = dplyr::filter(df_all, loglinear==TRUE, model!="adult", name!="Longspine thornyhead") %>%
  ggplot(aes(name,trend)) +
  geom_pointrange(aes(ymin=trend-2*trend_se, ymax=trend+2*trend_se),col="darkblue",size=0.8) +
  xlab("Species") + ylab("Trend in juvenile spatiotemporal sd (+/- 2SE)") +
  geom_hline(aes(yintercept=0),col="red",alpha=0.6) +
  coord_flip() +
  theme_bw()
#p2

# put both adults and juveniles on the same figure
p3 = dplyr::filter(df_all, loglinear==TRUE, name!="Longspine thornyhead") %>%
  dplyr::rename(Data=model) %>%
  ggplot(aes(name,trend,group=Data,col=Data)) +
  geom_pointrange(aes(ymin=trend-2*trend_se, ymax=trend+2*trend_se),size=0.8,alpha=0.8,
    position = position_dodge(width = 0.9)) +
  xlab("Species") +
  ylab("Trend in juvenile spatiotemporal sd (+/- 2SE)") +
  geom_hline(aes(yintercept=0),col="red",alpha=0.6) +
  coord_flip() +
  theme_bw() +
  scale_color_viridis(discrete=TRUE,end=0.8)
p3

dev.off()


aic_summary = dplyr::group_by(df_all, name, model) %>%
  dplyr::summarize(min_aic = min(aic,na.rm=T),
    diff_aic = aic[which(loglinear==TRUE)] - min_aic) %>%
  dplyr::arrange(model, name) %>%
  as.data.frame()
write.csv(aic_summary, "output/aic_summary.csv")



# Species of interest
species = read.csv("survey_data/species_list.csv")
names(species) = tolower(names(species))
species = dplyr::rename(species,
                        common_name = common.name,
                        scientific_name = scientific.name)
species$coef = NA
species$se = NA
species_df = species

for(i in 1:nrow(species_df)){
  comm_name = species_df$common_name[i]
  load(file=paste0("output/", sub(" ", "_", comm_name),"_all_models.RData"))
  indx = grep("b_epsilon_logit", names(ad_fit_ll$sd_report$value))
  species_df$coef[i] = ad_fit_ll$sd_report$value[indx]
  species_df$se[i] = ad_fit_ll$sd_report$sd[indx]
}

p4 = dplyr::filter(species_df,
              species %in% c("Greenspotted rockfish","Longspine thornyhead","Deepsea sole")==FALSE) %>%
  ggplot(aes(species,coef)) +
  geom_pointrange(aes(ymin=coef-2*se, ymax=coef+2*se),size=0.8,alpha=0.8,
                  position = position_dodge(width = 0.9)) +
  xlab("Species") +
  ylab("Trend in spatiotemporal sd (+/- 2SE)") +
  geom_hline(aes(yintercept=0),col="red",alpha=0.6) +
  coord_flip() +
  theme_bw() +
  scale_color_viridis(discrete=TRUE,end=0.8)




# Presence absence models
species = read.csv("survey_data/species_list.csv")
names(species) = tolower(names(species))
species = dplyr::rename(species,
                        common_name = common.name,
                        scientific_name = scientific.name)

for(i in 1:nrow(species)){
  print(i)
  comm_name = species$common_name[i]
  load(file=paste0("output/", sub(" ", "_", comm_name),"_presence_models.RData"))

  df = data.frame(name = comm_name,
                  model = c("adult", "adult"),
                  loglinear = c(FALSE,TRUE))
  df$aic = NA
  df$trend = NA
  df$trend_se = NA

  if(class(ad_fit)!="try-error") {
    df$aic[1] = AIC(ad_fit)
  }
  if(class(ad_fit_ll)!="try-error") {
    df$aic[2] = AIC(ad_fit_ll)
    df$trend[2] = ad_fit_ll$sd_report$value[which(names(ad_fit_ll$sd_report$value) == "b_epsilon_logit")]
    df$trend_se[2] = ad_fit_ll$sd_report$sd[which(names(ad_fit_ll$sd_report$value) == "b_epsilon_logit")]
  }
  if(i==1) {
    df_all = df
  } else {
    df_all = rbind(df_all, df)
  }
}

# plot trend coefficient for presence-absence models
p3 = dplyr::filter(df_all, loglinear==TRUE) %>%
  ggplot(aes(name,trend)) +
  geom_pointrange(aes(ymin=trend-2*trend_se, ymax=trend+2*trend_se),size=0.8,alpha=0.8,
                  position = position_dodge(width = 0.9)) +
  xlab("Species") +
  ylab("Trend in total presence-absence spatiotemporal sd (+/- 2SE)") +
  geom_hline(aes(yintercept=0),col="red",alpha=0.6) +
  ylim(-5,5) +
  coord_flip() +
  theme_bw() +
  scale_color_viridis(discrete=TRUE,end=0.8)

aic_summary = dplyr::group_by(df_all, name, model) %>%
  dplyr::summarize(min_aic = min(aic,na.rm=T),
                   diff_aic = aic[which(loglinear==TRUE)] - min_aic) %>%
  dplyr::arrange(model, name) %>%
  as.data.frame()
write.csv(aic_summary, "output/aic_summary_presence.csv")


# Positivie models
species = read.csv("survey_data/species_list.csv")
names(species) = tolower(names(species))
species = dplyr::rename(species,
                        common_name = common.name,
                        scientific_name = scientific.name)

for(i in 1:nrow(species)){
  print(i)
  comm_name = species$common_name[i]
  load(file=paste0("output/", sub(" ", "_", comm_name),"_positive_models.RData"))

  df = data.frame(name = comm_name,
                  model = c("adult", "adult"),
                  loglinear = c(FALSE,TRUE))
  df$aic = NA
  df$trend = NA
  df$trend_se = NA

  if(class(ad_fit)!="try-error") {
    df$aic[1] = AIC(ad_fit)
  }
  if(class(ad_fit_ll)!="try-error") {
    df$aic[2] = AIC(ad_fit_ll)
    df$trend[2] = ad_fit_ll$sd_report$value[which(names(ad_fit_ll$sd_report$value) == "b_epsilon_logit")]
    df$trend_se[2] = ad_fit_ll$sd_report$sd[which(names(ad_fit_ll$sd_report$value) == "b_epsilon_logit")]
  }
  if(i==1) {
    df_all = df
  } else {
    df_all = rbind(df_all, df)
  }
}

# plot trend coefficient for presence-absence models
p3 = dplyr::filter(df_all, loglinear==TRUE) %>%
  ggplot(aes(name,trend)) +
  geom_pointrange(aes(ymin=trend-2*trend_se, ymax=trend+2*trend_se),size=0.8,alpha=0.8,
                  position = position_dodge(width = 0.9)) +
  xlab("Species") +
  ylab("Trend in total positive spatiotemporal sd (+/- 2SE)") +
  geom_hline(aes(yintercept=0),col="red",alpha=0.6) +
  ylim(-5,5) +
  coord_flip() +
  theme_bw() +
  scale_color_viridis(discrete=TRUE,end=0.8)

aic_summary = dplyr::group_by(df_all, name, model) %>%
  dplyr::summarize(min_aic = min(aic,na.rm=T),
                   diff_aic = aic[which(loglinear==TRUE)] - min_aic) %>%
  dplyr::arrange(model, name) %>%
  as.data.frame()
write.csv(aic_summary, "output/aic_summary_positive.csv")
