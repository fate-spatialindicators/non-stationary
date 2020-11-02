library(sdmTMB)
library(sp)
library(ggplot2)
species = read.csv("survey_data/species_list.csv")
bio = readRDS("survey_data/wcbts_bio_2019-08-01.rds")
haul = readRDS("survey_data/wcbts_haul_2019-08-01.rds")
catch = readRDS("survey_data/wcbts_catch_2019-08-01.rds")
names(catch) = tolower(names(catch))

catch$trawl_id = as.character(catch$trawl_id)
haul$trawl_id = as.character(haul$trawl_id)
dat = dplyr::left_join(catch, haul)
# filter species from list
dat = dplyr::filter(dat, scientific_name %in% species$Scientific.Name)
# convert to UTM
coordinates(dat) <- c("longitude_dd", "latitude_dd")
proj4string(dat) <- CRS("+proj=longlat +datum=WGS84")
newproj = paste("+proj=utm +zone=10 ellps=WGS84")
dat <- spTransform(dat, CRS(newproj))
dat = as.data.frame(dat)
dat$lon = dat$longitude_dd/1000
dat$lat = dat$latitude_dd/1000
#dat$year = as.numeric(substr(dat$date_yyyymmdd,1,4))

models = list()
models_no_trend = list()
for(i in 1:length(unique(dat$scientific_name))) {
  
  sub = dplyr::filter(dat, scientific_name==unique(dat$scientific_name)[i])
  spde <- make_mesh(sub, c("lon", "lat"), cutoff = 20) 
  sub$depth_scaled = as.numeric(scale(sub$depth_m))
  sub$depth_scaled2 = sub$depth_scaled^2
  print(i)
  m = sdmTMB(cpue_kg_km2 ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
    data = sub, time = "year", spde = spde, family = tweedie(link = "log"),
    epsilon_model = "loglinear")
  models[[i]] = m
  
  m = sdmTMB(cpue_kg_km2 ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
    data = sub, time = "year", spde = spde, family = tweedie(link = "log"))
  models_no_trend[[i]] = m
}

# using ~ darblotched as an example (20)
par_trend = tidy(models[[20]])
par_trend = par_trend[grep("as.factor",par_trend$term),]
par_trend$year = 2003:2018
par_trend$model = "trend"
par_notrend = tidy(models_no_trend[[20]])
par_notrend = par_notrend[grep("as.factor",par_notrend$term),]
par_notrend$year = 2003:2018 + 0.2
par_notrend$model = "no_trend"
df = rbind(par_trend, par_notrend)
df$model =as.factor(df$model)
pdf("fig_3_example.pdf")
ggplot(df, aes(year, estimate,group=model,fill=model,col=model)) + 
  geom_pointrange(aes(ymin=estimate-std.error,ymax=estimate+std.error))
dev.off()

# summarize resulsts
out = data.frame("scientific_name" = unique(dat$scientific_name),stringsAsFactors = FALSE)
species_table = dplyr::group_by(bio, scientific_name) %>% 
  dplyr::summarise(common_name = common_name[1])
out = dplyr::left_join(out, species_table)
out$common_name = as.character(out$common_name)
out$b_est = NA
out$b_low = NA
out$b_hi = NA
out$sigma_E0 = NA
out$b_epsilon = NA
for(i in 1:nrow(out)) {
  m = models[[i]]
  out$b_epsilon[i] = m$tmb_obj$report()$b_epsilon
  out$sigma_E0[i] = m$tmb_obj$report()$sigma_E[1]
  out$b_est[i] = m$sd_report$par.fixed["b_epsilon_logit"]
  out$b_low[i] = m$sd_report$par.fixed["b_epsilon_logit"] - 
    1.96 * sqrt(m$sd_report$cov.fixed["b_epsilon_logit","b_epsilon_logit"])
  out$b_hi[i] = m$sd_report$par.fixed["b_epsilon_logit"] +
    1.96 * sqrt(m$sd_report$cov.fixed["b_epsilon_logit","b_epsilon_logit"])
}

out$b_est = 2 * exp(out$b_est) / (1+exp(out$b_est)) - 1
out$b_low = 2 * exp(out$b_low) / (1+exp(out$b_low)) - 1
out$b_hi = 2 * exp(out$b_hi) / (1+exp(out$b_hi)) - 1

pdf("fig_2_example_a.pdf")
ggplot(out, aes(common_name, b_est)) + 
  geom_pointrange(aes(ymin=b_low,ymax=b_hi)) + ylim(-0.5,0.5) + 
  coord_flip() + ylab("Slope on epsilon_st (log space)") + xlab("Common name")
dev.off()

df = expand.grid(common_name = out$common_name, year = 2003:2018)
df = left_join(df, out[,c("common_name","b_epsilon","sigma_E0")])
df$pred = exp(log(df$sigma_E0) + (df$year-2003+1) * df$b_epsilon)
df$common_name = as.factor(df$common_name)
pdf("fig_2_example_b.pdf")
ggplot(df, aes(year, pred,group=common_name,col=common_name)) + geom_line() + 
  xlab("Year") + ylab("Predicted epsilon(t)")
dev.off()

# perhaps do some kind of diagnostic, e.g. -- again darkblotched
i=20
sub_trend = dplyr::filter(dat, scientific_name==unique(dat$scientific_name)[i])
sub_trend$pred = predict(models[[i]])
sub_trend$model ="trend"
sub = dplyr::filter(dat, scientific_name==unique(dat$scientific_name)[i])
sub$pred = predict(models_no_trend[[i]])
sub$model ="no trend"
sub = dplyr::select(sub, year, cpue_kg_km2,model, pred)
sub_trend = dplyr::select(sub_trend, year, cpue_kg_km2,model, pred)

sub_bind = data.frame("year" = c(sub$year,sub_trend$year),
  "cpue_kg_km2" = c(sub$cpue_kg_km2,sub_trend$cpue_kg_km2),
  "model" = c(sub$model, sub_trend$model),
  pred = c(sub$pred$est,sub_trend$pred$est))
pdf("fig_s1.pdf")
dplyr::group_by(dplyr::filter(sub_bind, cpue_kg_km2<200,cpue_kg_km2>0), year, model) %>% 
  dplyr::summarise(rho = cor(pred, log(cpue_kg_km2))) %>% 
  as.data.frame()
dev.off()

pdf("fig_s2.pdf")
ggplot(dplyr::filter(sub_bind, cpue_kg_km2>0), 
  aes(pred, log(cpue_kg_km2),col=model,group=model)) + 
  geom_point(alpha=0.3,size=0.4) + 
  facet_wrap(~year) + xlab("Predictions") + 
  ylab("Log(observed CPUE)")
dev.off()

i=20
sub = dplyr::filter(dat, scientific_name==unique(dat$scientific_name)[i])
sub$pred_trend = predict(models[[i]])$est
sub$pred = predict(models_no_trend[[i]])$est
ggplot(dplyr::filter(sub, cpue_kg_km2>0), 
  aes(pred, pred_trend)) + 
  geom_point(alpha=0.3,size=0.4) + 
  facet_wrap(~year) + xlab("Predictions - no trend") + 
  ylab("Predictions - trend")

# fit same models to temp and o2
haul = readRDS("survey_data/wcbts_haul_2019-08-01.rds")

coordinates(haul) <- c("longitude_dd", "latitude_dd")
proj4string(haul) <- CRS("+proj=longlat +datum=WGS84")
newproj = paste("+proj=utm +zone=10 ellps=WGS84")
haul <- spTransform(haul, CRS(newproj))
haul = as.data.frame(haul)
haul$lon = haul$longitude_dd/1000
haul$lat = haul$latitude_dd/1000
haul$depth_scaled = as.numeric(scale(haul$depth_hi_prec_m))
haul$depth_scaled2 = haul$depth_scaled^2
haul$year = as.numeric(substr(haul$date,1,4))
haul$month = as.numeric(substr(haul$date,5,6))
haul$day = as.numeric(substr(haul$date,7,8))
haul$yday = mdy.date(haul$month, haul$day, haul$year)- mdy.date(1, 1, haul$year)

# First fit model using observed temp as response
spde <- make_mesh(haul, c("lon", "lat"), cutoff = 20) 
fit_temp = sdmTMB(temperature_at_gear_c_der ~ 0 + depth_scaled + depth_scaled2 + 
    as.factor(year) + yday + I(yday^2),
  data = haul, time = "year", spde = spde, epsilon_model = "loglinear")
fit_temp$sd_report 

# Second fit model using observed o2 as response, subset of years
o2_haul = dplyr::filter(haul, year%in%seq(2010,2015))
spde <- make_mesh(o2_haul, c("lon", "lat"), cutoff = 20) 
fit_o2 = sdmTMB(o2_at_gear_ml_per_l_der ~ 0 + depth_scaled + depth_scaled2 + 
    as.factor(year) + yday + I(yday^2),
  data = o2_haul, 
  time = "year", spde = spde, epsilon_model = "loglinear")
fit_o2$sd_report
