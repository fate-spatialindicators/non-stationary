## Fit Stationary and Non-Stationary Models to zero and positive responses ##
## In parallel ##

remotes::install_github("pbs-assess/sdmTMB", "nonstationary") # v. 0.0.19.9001
library(sdmTMB)
library(dplyr)
library(future)
options(future.globals.maxSize = 4e3 * 1024 ^ 2) # ~4GB
plan(multisession, workers = 8) # or just specify # cores manually
if (!dir.exists("output")) dir.create("output")

# 15 -> ~ 600 knots; 20 -> 389 knots; 25 -> 294 knots; 30 -> 221 knots
n_cutoff <- 20
species <- read.csv("survey_data/species_list.csv", fileEncoding="UTF-8-BOM")
names(species) <- tolower(names(species))
species <- dplyr::rename(species,
                         common_name = species,
                         scientific_name = scientific.name
)

dat <- readRDS("survey_data/all_data.rds")

dat <- left_join(select(species, scientific_name, common_name), dat)

grid <- readRDS("data/wc_grid.rds")
grid <- rename(grid, lon = X, lat = Y)
grid$depth_scaled <- as.numeric(scale(grid$depth))
grid$depth_scaled2 <- grid$depth_scaled^2

dat$year <- as.factor(dat$year)
temp <- (dat$depth_m - mean(-grid$depth))
dat$depth_scaled <- temp / sd(grid$depth)
dat$depth_scaled2 <- dat$depth_scaled^2
dat$time <- as.numeric(as.character(dat$year)) - floor(mean(unique(as.numeric(as.character(dat$year)))))

#dat$time <- dat$time - mean(dat$time)
# # spatial + spatiotemporal:
# fit_models <- function(sub) {
#   spde <- make_mesh(sub, c("lon", "lat"), cutoff = n_cutoff, type = "cutoff")
#   ad_fit <- tryCatch({
#     sdmTMB(cpue_kg_km2 ~ 0 + depth_scaled + depth_scaled2 + year,
#            data = sub, time = "year", spde = spde, family = tweedie(link = "log"),
#            control = sdmTMBcontrol(nlminb_loops = 2, newton_loops = 1)
#     )}, error = function(e) NA)
#   #ad_fit <- refit_model_if_needed(ad_fit)
#   ad_fit_ll <- tryCatch({
#     sdmTMB(cpue_kg_km2 ~ 0 + depth_scaled + depth_scaled2 + year,
#            data = sub, time = "year", spde = spde,
#            family = tweedie(link = "log"), epsilon_predictor = "time",
#            control = sdmTMBcontrol(nlminb_loops = 2, newton_loops = 1)
#     )}, error = function(e) NA)
#   #ad_fit_ll <- refit_model_if_needed(ad_fit_ll)
#   save(ad_fit, ad_fit_ll,
#        file = paste0("output/", sub(" ", "_", sub$common_name[[1]]), "_all_models.RData")
#   )
# }

#split(dat, dat$scientific_name) %>%
#  furrr::future_walk(fit_models)
# purrr::walk(fit_models) # serial for testing

# AR1 spatiotemporal:
fit_models_ar1 <- function(sub) {
  spde <- make_mesh(sub, c("lon", "lat"), cutoff = n_cutoff, type = "cutoff")
  sub$fyear <- as.factor(sub$year)
  ad_fit <- tryCatch({
    sdmTMB(cpue_kg_km2 ~ 0 + depth_scaled + depth_scaled2 + fyear,
           spatiotemporal = "ar1",
           spatial = "off",
           data = sub,
           time = "year",
           mesh = spde,
           family = tweedie(link = "log"),
           control = sdmTMBcontrol(nlminb_loops = 2, newton_loops = 1),
           priors = sdmTMBpriors(
             #phi = halfnormal(0, 10),
             #tweedie_p = normal(1.5, 2),
             #ar1_rho = normal(0.5, 0.1),
             matern_st = pc_matern(range_gt = 10, sigma_lt = 5))
    )}, error = function(e) NA)
  #ad_fit <- refit_model_if_needed(ad_fit)
  ad_fit_ll <- tryCatch({
    sdmTMB(cpue_kg_km2 ~ 0 + depth_scaled + depth_scaled2 + fyear,
           spatiotemporal = "ar1",
           spatial = "off",
           data = sub,
           time = "year",
           mesh = spde,
           family = tweedie(link = "log"),
           experimental = list(epsilon_predictor = "time", epsilon_model = "trend"),
           control = sdmTMBcontrol(nlminb_loops = 2, newton_loops = 1,
                                   lower = list(b_epsilon = -1),
                                   upper = list(b_epsilon = 1)),
           priors = sdmTMBpriors(
             #phi = halfnormal(0, 10),
             #tweedie_p = normal(1.5, 2),
             #ar1_rho = normal(0.5, 0.1),
             matern_st = pc_matern(range_gt = 10, sigma_lt = 5))
    )}, error = function(e) NA)
  #ad_fit_ll <- refit_model_if_needed(ad_fit_ll)
  save(ad_fit, ad_fit_ll,
       file = paste0("output/", sub(" ", "_", sub$common_name[[1]]), "_ar1_priors.RData")
  )
}

# refit_model_if_needed <- function(m) {
#   if (class(m)!="try-error") {
#     if (max(m$gradients) > 0.01) {
#       m <- tryCatch({
#         sdmTMB::run_extra_optimization(m,
#                                        nlminb_loops = 1L,
#                                        newton_steps = 1L
#         )
#       }, error = function(e) m)
#     }
#   }
#   m
# }

df = dplyr::group_by(dat, scientific_name) %>%
  dplyr::summarise(common_name = common_name[1])

species_names = saveRDS(df, file="species.rds")
split(dat, dat$scientific_name) %>%
  furrr::future_walk(fit_models_ar1)
# purrr::walk(fit_models_ar1) # serial for testing

plan(sequential) # avoid RStudio crashes on Session Restart R
