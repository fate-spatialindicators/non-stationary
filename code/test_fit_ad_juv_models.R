devtools::install_github("pbs-assess/sdmTMB","epsilon-heterogeneity")
library(sdmTMB)
library(sp)
library(ggplot2)
library(ggsidekick)

# make plots of variability in residuals for models of positive catch rates
dat = readRDS("Sebastes crameri_juv_cpue.rds")

sub = dplyr::filter(dat, juv_cpue_kg_km2 > 0)
sub$depth_scaled = as.numeric(scale(sub$depth_m))
sub$depth_scaled2 = sub$depth_scaled^2
spde <- make_mesh(sub, c("lon", "lat"), cutoff = 10)
# check knots with spde$mesh$n
pos_juv_fit = sdmTMB(juv_cpue_kg_km2 ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
  data = sub, time = "year", spde = spde, family = Gamma(link = "log"))
sub$resid = sub$juv_cpue_kg_km2 - predict(pos_juv_fit)$est
# plot of residual variation over time
g1 = dplyr::group_by(sub, year) %>%
  dplyr::summarize(sd = sd(resid)) %>%
  ggplot(aes(year, log(sd))) + geom_point() + ylab("Ln(SD) residuals juv cpue")

sub = dplyr::filter(dat, adult_cpue_kg_km2 > 0)
sub$depth_scaled = as.numeric(scale(sub$depth_m))
sub$depth_scaled2 = sub$depth_scaled^2
spde <- make_mesh(sub, c("lon", "lat"), cutoff = 10)
pos_ad_fit = sdmTMB(adult_cpue_kg_km2 ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
  data = sub, time = "year", spde = spde, family = Gamma(link = "log"))
sub$resid = sub$adult_cpue_kg_km2 - predict(pos_ad_fit)$est
# plot of residual variation over time
g2 = dplyr::group_by(sub, year) %>%
  dplyr::summarize(sd = sd(resid)) %>%
  ggplot(aes(year, log(sd))) + geom_point() + ylab("Ln(SD) residuals ad cpue")

# fit tweedie model to total cpue
sub$depth_scaled = as.numeric(scale(sub$depth_m))
sub$depth_scaled2 = sub$depth_scaled^2
spde <- make_mesh(sub, c("lon", "lat"), cutoff = 10)
juv_tweedie = sdmTMB(juv_cpue_kg_km2 ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
  data = sub, time = "year", spde = spde, family = tweedie(link = "log"),
  epsilon_model = "loglinear")
