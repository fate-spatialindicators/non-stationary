## Fit Stationary and Non-Stationary Models to zero and positive responses ##

# new data,
# smaller knots, 600,
# delta lognormal instead of gamma
#remotes::install_github("pbs-assess/sdmTMB", ref = "epsilon-heterogeneity")
library(sdmTMB)
library(dplyr)

n_cutoff = 15 # generates ~ 600 knots
# Species of interest
species = read.csv("survey_data/species_list.csv")
names(species) = tolower(names(species))
species = dplyr::rename(species,
                        common_name = common.name,
                        scientific_name = scientific.name)

dat = readRDS("survey_data/all_data.rds")
# loop over species to fit models ----
for(i in 1:nrow(species)){

  comm_name = species$common_name[i]

  sub = dplyr::filter(dat, scientific_name == species$scientific_name[i])

  sub$depth_scaled = as.numeric(scale(sub$depth_m))
  sub$depth_scaled2 = sub$depth_scaled^2

  spde = make_mesh(sub, c("lon", "lat"), cutoff = n_cutoff)

  # total cpue
  ad_fit = try(sdmTMB(cpue_kg_km2 ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
                   data = sub, time = "year", spde = spde, family = tweedie(link = "log"))
               )

  # scale time
  sub$time = sub$year - min(sub$year) + 1
  ad_fit_ll = try(sdmTMB(cpue_kg_km2 ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
                      data = sub, time = "year", spde = spde,
                      family = tweedie(link = "log"), epsilon_predictor = "time")
  )

  # calculate avg bottom temp
  temp = dplyr::group_by(sub, year) %>%
    dplyr::summarize(mean_temp = mean(temperature_at_gear_c_der,na.rm=T))
  temp$mean_temp = scale(temp$mean_temp)
  sub = dplyr::left_join(sub, temp)
  ad_fit_ll_temp = try(sdmTMB(cpue_kg_km2 ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
                         data = sub, time = "year", spde = spde,
                         family = tweedie(link = "log"), epsilon_predictor = "mean_temp")
  )

  # drop out high memory objects
  ad_fit$data = NULL
  ad_fit$tmb_data = NULL
  ad_fit$spde = NULL
  ad_fit_ll$data = NULL
  ad_fit_ll$tmb_data = NULL
  ad_fit_ll$spde = NULL
  ad_fit_ll_temp$data = NULL
  ad_fit_ll_temp$tmb_data = NULL
  ad_fit_ll_temp$spde = NULL
  if (!dir.exists("output")) {dir.create("output")}
  save(ad_fit, ad_fit_ll, ad_fit_ll_temp,file=paste0("output/", sub(" ", "_", comm_name),"_all_models.RData"))
}

# fit presence - absence models
for(i in 1:nrow(species)){

  comm_name = species$common_name[i]

  sub = dplyr::filter(dat, scientific_name == species$scientific_name[i])

  sub$depth_scaled = as.numeric(scale(sub$depth_m))
  sub$depth_scaled2 = sub$depth_scaled^2

  spde = make_mesh(sub, c("lon", "lat"), cutoff = n_cutoff)

  sub$presence = ifelse(sub$cpue_kg_km2 > 0, 1, 0)
  # total cpue
  ad_fit = try(sdmTMB(presence ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
                      data = sub, time = "year", spde = spde, family = binomial(link = "logit"))
  )

  # scale time
  sub$time = sub$year - min(sub$year) + 1
  ad_fit_ll = try(sdmTMB(presence ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
                         data = sub, time = "year", spde = spde, family = binomial(link = "logit"),
                         epsilon_predictor = "time")
  )

  # calculate avg bottom temp
  temp = dplyr::group_by(sub, year) %>%
    dplyr::summarize(mean_temp = mean(temperature_at_gear_c_der,na.rm=T))
  temp$mean_temp = scale(temp$mean_temp)
  sub = dplyr::left_join(sub, temp)
  ad_fit_ll_temp = try(sdmTMB(presence ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
                              data = sub, time = "year", spde = spde,
                              family = binomial(link = "logit"), epsilon_predictor = "mean_temp")
  )

  # drop out high memory objects
  ad_fit$data = NULL
  ad_fit$tmb_data = NULL
  ad_fit$spde = NULL
  ad_fit_ll$data = NULL
  ad_fit_ll$tmb_data = NULL
  ad_fit_ll$spde = NULL
  ad_fit_ll_temp$data = NULL
  ad_fit_ll_temp$tmb_data = NULL
  ad_fit_ll_temp$spde = NULL
  if (!dir.exists("output")) {dir.create("output")}
  save(ad_fit, ad_fit_ll, ad_fit_ll_temp,file=paste0("output/", sub(" ", "_", comm_name),"_presence_models.RData"))
}


# fit positive models
for(i in 1:nrow(species)){

  comm_name = species$common_name[i]

  sub = dplyr::filter(dat, scientific_name == species$scientific_name[i])

  sub$depth_scaled = as.numeric(scale(sub$depth_m))
  sub$depth_scaled2 = sub$depth_scaled^2

  sub = dplyr::filter(sub, cpue_kg_km2 > 0)
  spde = make_mesh(sub, c("lon", "lat"), cutoff = n_cutoff)

  #sub$presence = ifelse(sub$cpue_kg_km2 > 0, 1, 0)
  # total cpue
  ad_fit = try(sdmTMB(cpue_kg_km2 ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
                      data = sub, time = "year", spde = spde, family = lognormal(link = "log"))
  )

  sub$time = sub$year - min(sub$year) + 1
  ad_fit_ll = try(sdmTMB(cpue_kg_km2 ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
                         data = sub, time = "year", spde = spde, family = lognormal(link = "log"),
                         epsilon_predictor = "time")
  )

  # calculate avg bottom temp
  temp = dplyr::group_by(sub, year) %>%
    dplyr::summarize(mean_temp = mean(temperature_at_gear_c_der,na.rm=T))
  temp$mean_temp = scale(temp$mean_temp)
  sub = dplyr::left_join(sub, temp)
  ad_fit_ll_temp = try(sdmTMB(presence ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
                              data = sub, time = "year", spde = spde,
                              family = lognormal(link = "log"), epsilon_predictor = "mean_temp")
  )

  # drop out high memory objects
  ad_fit$data = NULL
  ad_fit$tmb_data = NULL
  ad_fit$spde = NULL
  ad_fit_ll$data = NULL
  ad_fit_ll$tmb_data = NULL
  ad_fit_ll$spde = NULL
  ad_fit_ll_temp$data = NULL
  ad_fit_ll_temp$tmb_data = NULL
  ad_fit_ll_temp$spde = NULL
  if (!dir.exists("output")) {dir.create("output")}
  save(ad_fit, ad_fit_ll, ad_fit_ll_temp,file=paste0("output/", sub(" ", "_", comm_name),"_positive_models.RData"))
}
