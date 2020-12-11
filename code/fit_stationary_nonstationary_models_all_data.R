## Fit Stationary and Non-Stationary Models to zero and positive responses ##

#remotes::install_github("pbs-assess/sdmTMB", ref = "epsilon-heterogeneity")
library(sdmTMB)
library(dplyr)

# Species of interest
species = read.csv("survey_data/species_list_revised.csv")
names(species) = tolower(names(species))
species = dplyr::rename(species,
                        common_name = common.name,
                        scientific_name = scientific.name,
                        juv_threshold = max.length.cm)

# loop over species to fit models ----
for(i in 1:nrow(species)){

  comm_name = species$common_name[i]
  dat = readRDS(paste0("data/", sub(" ", "_", comm_name), "_expanded.rds"))

  dat$depth_scaled = as.numeric(scale(dat$depth_m))
  dat$depth_scaled2 = dat$depth_scaled^2

  spde = make_mesh(dat, c("lon", "lat"), cutoff = 10)

  # adults and juveniles combined
  total_fit = sdmTMB(cpue_kg_km2 ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
                       data = dat, time = "year", spde = spde, family = tweedie(link = "log"),
                       epsilon_model = "loglinear")
  # juveniles
  juv_fit = sdmTMB(juv_cpue_kg_km2 ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
                     data = dat, time = "year", spde = spde, family = tweedie(link = "log"),
                     epsilon_model = "loglinear")
  # adults
  ad_fit = sdmTMB(adult_cpue_kg_km2 ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
                   data = dat, time = "year", spde = spde, family = tweedie(link = "log"),
                   epsilon_model = "loglinear")

  # same as the above 3 models for the stationary case....
  total_fit_stationary = sdmTMB(cpue_kg_km2 ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
                     data = dat, time = "year", spde = spde, family = tweedie(link = "log"))
  juv_fit_stationary = sdmTMB(juv_cpue_kg_km2 ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
                   data = dat, time = "year", spde = spde, family = tweedie(link = "log"))
  ad_fit_stationary = sdmTMB(adult_cpue_kg_km2 ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
                  data = dat, time = "year", spde = spde, family = tweedie(link = "log"))

  rm(dat)
  if (!dir.exists("output")) {dir.create("output")}
  save.image(file=paste0("output/", sub(" ", "_", comm_name),"_all_models.RData"))
  print(paste0("species ", i, " of ", nrow(species)), "completed")
}
