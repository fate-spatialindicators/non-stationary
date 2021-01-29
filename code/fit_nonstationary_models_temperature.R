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

# filter out a single species -- all we need is temp data
sub = dplyr::filter(dat,
                    scientific_name == unique(dat$scientific_name)[1])

sub$depth_scaled = as.numeric(scale(sub$depth_m))
sub$depth_scaled2 = sub$depth_scaled^2

sub = dplyr::filter(sub, !is.na(temperature_at_gear_c_der))
spde = make_mesh(sub, c("lon", "lat"), cutoff = n_cutoff)

# total cpue
fit = try(sdmTMB(temperature_at_gear_c_der ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
                    data = sub, time = "year", spde = spde)
)

# scale time
sub$time = sub$year - min(sub$year) + 1
ad_fit_ll = try(sdmTMB(cpue_kg_km2 ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
                       data = sub, time = "year", spde = spde, epsilon_predictor = "time")
)

# drop out high memory objects
ad_fit$data = NULL
ad_fit$tmb_data = NULL
ad_fit$spde = NULL
ad_fit_ll$data = NULL
ad_fit_ll$tmb_data = NULL
ad_fit_ll$spde = NULL

if (!dir.exists("output")) {dir.create("output")}
save.image(file=paste0("output/temp_all_models.RData"))
