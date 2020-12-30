## Fit Stationary and Non-Stationary Models to zero and positive responses ##

remotes::install_github("pbs-assess/sdmTMB", ref = "epsilon-heterogeneity")
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

  spde = make_mesh(dat, c("lon", "lat"), cutoff = 25) # 20 -> 389 knots; 25 -> 294 knots; 30 -> 221 knots

  # stationary fits
  m_total = try(sdmTMB(cpue_kg_km2 ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
                       data = dat, time = "year", spde = spde, family = tweedie(link = "log"),
                       nlminb_loops = 2, newton_steps = 1)
  )
  m_juv = try(sdmTMB(juv_cpue_kg_km2 ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
                     data = dat, time = "year", spde = spde, family = tweedie(link = "log"),
                     nlminb_loops = 2, newton_steps = 1)
  )
  m_adult = try(sdmTMB(adult_cpue_kg_km2 ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
                       data = dat, time = "year", spde = spde, family = tweedie(link = "log"),
                       nlminb_loops = 2, newton_steps = 1)
  )

  print(paste0("species ", i, " of ", nrow(species), " stationary complete"))

  # with loglinear on sigma_epsilon
  m_total_ll = try(sdmTMB(cpue_kg_km2 ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
                          data = dat, time = "year", spde = spde, family = tweedie(link = "log"),
                          nlminb_loops = 2, newton_steps = 1,
                          epsilon_model = "loglinear")
  )
  m_juv_ll = try(sdmTMB(juv_cpue_kg_km2 ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
                        data = dat, time = "year", spde = spde, family = tweedie(link = "log"),
                        nlminb_loops = 2, newton_steps = 1,
                        epsilon_model = "loglinear")
  )
  m_adult_ll = try(sdmTMB(adult_cpue_kg_km2 ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
                          data = dat, time = "year", spde = spde, family = tweedie(link = "log"),
                          nlminb_loops = 2, newton_steps = 1,
                          epsilon_model = "loglinear")
  )

  print(paste0("species ", i, " of ", nrow(species), " loglinear complete"))

  # # with random effects on sigma_epsilon
  # m_total_re = try(sdmTMB(cpue_kg_km2 ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
  #                         data = dat, time = "year", spde = spde, family = tweedie(link = "log"),
  #                         nlminb_loops = 2, newton_steps = 1,
  #                         epsilon_model = "re")
  # )
  # m_juv_re = try(sdmTMB(juv_cpue_kg_km2 ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
  #                       data = dat, time = "year", spde = spde, family = tweedie(link = "log"),
  #                       nlminb_loops = 2, newton_steps = 1,
  #                       epsilon_model = "re")
  # )
  # m_adult_re = try(sdmTMB(adult_cpue_kg_km2 ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
  #                         data = dat, time = "year", spde = spde, family = tweedie(link = "log"),
  #                         nlminb_loops = 2, newton_steps = 1,
  #                         epsilon_model = "re")
  # )
  #
  # print(paste0("species ", i, " of ", nrow(species), " re complete"))

  # # with combination of loglinear and random effects on sigma_epsilon
  # m_total_ll_re = try(sdmTMB(cpue_kg_km2 ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
  #                        data = dat, time = "year", spde = spde, family = tweedie(link = "log"),
  #                        nlminb_loops = 2, newton_steps = 1,
  #                        epsilon_model = "loglinear-re")
  # )
  # m_juv_ll_re = try(sdmTMB(juv_cpue_kg_km2 ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
  #                      data = dat, time = "year", spde = spde, family = tweedie(link = "log"),
  #                      nlminb_loops = 2, newton_steps = 1,
  #                      epsilon_model = "loglinear-re")
  # )
  # m_adult_ll_re = try(sdmTMB(adult_cpue_kg_km2 ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
  #                     data = dat, time = "year", spde = spde, family = tweedie(link = "log"),
  #                     nlminb_loops = 2, newton_steps = 1,
  #                     epsilon_model = "loglinear-re")
  # )
  #
  # print(paste0("species ", i, " of ", nrow(species), " loglinear/re complete"))

  rm(dat, spde)

  if (!dir.exists("output")) {dir.create("output")}
  save.image(file=paste0("output/", sub(" ", "_", comm_name),"_all_models.RData"))

  try(rm(m_adult_re,m_juv_re,m_total_re,m_adult_ll,m_juv_ll,m_total_ll,m_adult,m_juv,m_total), silent = TRUE)

}
