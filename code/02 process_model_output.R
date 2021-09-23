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

  load(file=paste0("output/", sub(" ", "_", comm_name),"_ar1_priors.RData"))

  df = data.frame(name = comm_name,
    model = c("adult", "adult"),
    loglinear = c(FALSE,TRUE))
  df$pred_dens = NA
  df$trend = NA
  df$trend_se = NA

  if(class(ad_fit)!="try-error") {
    df$pred_dens[1] = ad_fit$sum_loglik
  }
  if(class(ad_fit_ll)!="try-error") {
    df$pred_dens[2] = ad_fit_ll$sum_loglik
    df_trend = data.frame(trend = rep(NA, length(ad_fit_ll$models)),
                          trend_se = rep(NA, length(ad_fit_ll$models)))
    for(k in 1:length(ad_fit_ll$models)) {
      df_trend$trend[k] = ad_fit_ll$models[[k]]$sd_report$value[which(names(ad_fit_ll$models[[k]]$sd_report$value) == "b_epsilon")]
      df_trend$trend_se[k] = ad_fit_ll$models[[k]]$sd_report$sd[which(names(ad_fit_ll$models[[k]]$sd_report$value) == "b_epsilon")]
    }
    df_trend$inv_var = 1/df_trend$trend_se^2
    df$trend[2] = sum(df_trend$trend * df_trend$inv_var) / sum(df_trend$inv_var)
    df$trend_se[2] = sqrt(1/sum(df_trend$inv_var))
  }

  if(i==1) {
    df_all = df
  } else {
    df_all = rbind(df_all, df)
  }
}

# summarize the predictive density diffs
dens_summary = dplyr::group_by(df_all, name) %>%
  dplyr::summarize(max_diff = max(pred_dens,na.rm=T),
                   diff_dens = max_diff - pred_dens[which(loglinear==TRUE)]) %>%
  dplyr::arrange(name) %>%
  as.data.frame()
write.csv(dens_summary, "output/pred_dens_summary.csv")

df_all = left_join(df_all, dens_summary)
df_all$model = ifelse(df_all$diff_dens == 0, "trend", "no trend")
# look at coefficients for adult models
jpeg("plots/Trend_summaries.jpeg")

p1 = dplyr::filter(df_all, loglinear==TRUE) %>%
  ggplot(aes(name,trend, col=model)) +
  geom_pointrange(aes(ymin=trend-2*trend_se, ymax=trend+2*trend_se),size=0.8,alpha=0.6) +
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
  geom_smooth(method="lm")# +
  #xlim(0.25,1) +
  #ylim(-0.25,0.11)
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






