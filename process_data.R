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
    df$trend[4] = m_adult_ll$sd_report$value[which(names(m_adult_ll$sd_report$value) == "b_epsilon_logit")]
    df$trend_se[4] = m_adult_ll$sd_report$sd[which(names(m_adult_ll$sd_report$value) == "b_epsilon_logit")]
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
  geom_pointrange(aes(ymin=trend-2*trend_se, ymax=trend+2*trend_se)) +
  xlab("Species") + ylab("Trend (+/- 2SE)") +
  geom_hline(aes(yintercept=0),col="red",alpha=0.6) +
  coord_flip() +
  ggtitle("Trends in adult spatiotemporal sd")
p1

p2 = dplyr::filter(df_all, loglinear==TRUE, model!="adult", name!="Longspine thornyhead") %>%
  ggplot(aes(name,trend)) +
  geom_pointrange(aes(ymin=trend-2*trend_se, ymax=trend+2*trend_se)) +
  xlab("Species") + ylab("Trend (+/- 2SE)") +
  geom_hline(aes(yintercept=0),col="red",alpha=0.6) +
  coord_flip() +
  ggtitle("Trends in juvenile spatiotemporal sd")
p2

dev.off()


aic_summary = dplyr::group_by(df_all, name, model) %>%
  dplyr::summarize(min_aic = min(aic,na.rm=T),
    diff_aic = aic[which(loglinear==TRUE)] - min_aic) %>%
  dplyr::arrange(model, name) %>%
  as.data.frame()
write.csv(aic_summary, "output/aic_summary.csv")
