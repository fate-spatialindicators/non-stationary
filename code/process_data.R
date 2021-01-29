library(sdmTMB)
library(dplyr)
library(ggplot2)
library(viridis)

# Species of interest
species = read.csv("survey_data/species_list_revised.csv")
names(species) = tolower(names(species))
species = dplyr::rename(species,
  common_name = common.name,
  scientific_name = scientific.name,
  juv_threshold = max.length.cm)

for(i in 1:nrow(species)){

  comm_name = species$common_name[i]
  dat = readRDS(paste0("data/", sub(" ", "_", comm_name), "_expanded.rds"))

  dat$depth_scaled = as.numeric(scale(dat$depth_m))
  dat$depth_scaled2 = dat$depth_scaled^2

  load(file=paste0("output/", sub(" ", "_", comm_name),"_all_models.RData"))

  df = data.frame(name = comm_name,
    model = c("adult", "adult","juvenile","juvenile"),
    loglinear = c(FALSE,TRUE,FALSE,TRUE))
  df$aic = NA
  df$trend = NA
  df$trend_se = NA

  if(class(m_adult)!="try-error") {
    df$aic[1] = AIC(m_adult)
  }
  if(class(m_adult_ll)!="try-error") {
    df$aic[2] = AIC(m_adult_ll)
    df$trend[2] = m_adult_ll$sd_report$value[which(names(m_adult_ll$sd_report$value) == "b_epsilon_logit")]
    df$trend_se[2] = m_adult_ll$sd_report$sd[which(names(m_adult_ll$sd_report$value) == "b_epsilon_logit")]
  }
  if(class(m_juv)!="try-error") {
    df$aic[3] = AIC(m_juv)
  }
  if(class(m_juv_ll)!="try-error") {
    df$aic[4] = AIC(m_juv_ll)
    df$trend[4] = m_juv_ll$sd_report$value[which(names(m_juv_ll$sd_report$value) == "b_epsilon_logit")]
    df$trend_se[4] = m_juv_ll$sd_report$sd[which(names(m_juv_ll$sd_report$value) == "b_epsilon_logit")]
  }

  if(i==1) {
    df_all = df
  } else {
    df_all = rbind(df_all, df)
  }
}


# look at coefficients for adult models
pdf("plots/Trend_summaries.pdf")

p1 = dplyr::filter(df_all, loglinear==TRUE, model=="adult", name!="Longspine thornyhead") %>%
  ggplot(aes(name,trend)) +
  geom_pointrange(aes(ymin=trend-2*trend_se, ymax=trend+2*trend_se),col="darkblue",size=0.8) +
  xlab("Species") + ylab("Trend in adult spatiotemporal sd (+/- 2SE)") +
  geom_hline(aes(yintercept=0),col="red",alpha=0.6) +
  coord_flip() +
  theme_bw()
#p1

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
