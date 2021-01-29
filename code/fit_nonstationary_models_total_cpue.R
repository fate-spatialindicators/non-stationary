## Fit Stationary and Non-Stationary Models to zero and positive responses ##

#remotes::install_github("pbs-assess/sdmTMB", ref = "epsilon-heterogeneity")
library(sdmTMB)
library(dplyr)

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

  spde = make_mesh(sub, c("lon", "lat"), cutoff = 10)

  # total cpue
  ad_fit = try(sdmTMB(cpue_kg_km2 ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
                   data = sub, time = "year", spde = spde, family = tweedie(link = "log"),
                   nlminb_loops = 2, newton_steps = 1)
               )

  sub$time = sub$year - min(sub$year) + 1
  ad_fit_ll = try(sdmTMB(cpue_kg_km2 ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
                      data = sub, time = "year", spde = spde, family = tweedie(link = "log"),
                      nlminb_loops = 2, newton_steps = 1, epsilon_predictor = "time")
  )

  if (!dir.exists("output")) {dir.create("output")}
  save.image(file=paste0("output/", sub(" ", "_", comm_name),"_all_models.RData"))
}

# fit presence - absence models
for(i in 1:nrow(species)){

  comm_name = species$common_name[i]

  sub = dplyr::filter(dat, scientific_name == species$scientific_name[i])

  sub$depth_scaled = as.numeric(scale(sub$depth_m))
  sub$depth_scaled2 = sub$depth_scaled^2

  spde = make_mesh(sub, c("lon", "lat"), cutoff = 10)

  sub$presence = ifelse(sub$cpue_kg_km2 > 0, 1, 0)
  # total cpue
  ad_fit = try(sdmTMB(presence ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
                      data = sub, time = "year", spde = spde, family = binomial(link = "logit"))
  )

  sub$time = sub$year - min(sub$year) + 1
  ad_fit_ll = try(sdmTMB(presence ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
                         data = sub, time = "year", spde = spde, family = binomial(link = "logit"),
                         epsilon_predictor = "time")
  )

  if (!dir.exists("output")) {dir.create("output")}
  save.image(file=paste0("output/", sub(" ", "_", comm_name),"_presence_models.RData"))
}


# fit positive models
for(i in 1:nrow(species)){

  comm_name = species$common_name[i]

  sub = dplyr::filter(dat, scientific_name == species$scientific_name[i])

  sub$depth_scaled = as.numeric(scale(sub$depth_m))
  sub$depth_scaled2 = sub$depth_scaled^2

  sub = dplyr::filter(sub, cpue_kg_km2 > 0)
  spde = make_mesh(sub, c("lon", "lat"), cutoff = 10)

  #sub$presence = ifelse(sub$cpue_kg_km2 > 0, 1, 0)
  # total cpue
  ad_fit = try(sdmTMB(cpue_kg_km2 ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
                      data = sub, time = "year", spde = spde, family = Gamma(link = "log"))
  )

  sub$time = sub$year - min(sub$year) + 1
  ad_fit_ll = try(sdmTMB(cpue_kg_km2 ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
                         data = sub, time = "year", spde = spde, family = Gamma(link = "log"),
                         epsilon_predictor = "time")
  )

  if (!dir.exists("output")) {dir.create("output")}
  save.image(file=paste0("output/", sub(" ", "_", comm_name),"_positive_models.RData"))
}
