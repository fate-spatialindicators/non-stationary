library(sdmTMB)

# simulation will have baseline spatiotemporal variation of 0.2. So trend
# is interpreted relative to that
grid = expand.grid("rho" = c(0),
                   "b_trend" = c(0.01, 0.02, 0.04, 0.08),
                   "obs_sd" = c(0.1, 0.2, 0.4),
                   "iter" = 1:50)

n_site = 50
set.seed(42)
x <- runif(n_site, -1, 1)
y <- runif(n_site, -1, 1)
time_steps <- 10

df = expand.grid("site"=1:length(x), "time" = as.factor(1:time_steps))
X <- model.matrix(~ -1 + time, df)

loc <- data.frame(x = x, y = y)
mesh <- make_mesh(loc, xy_cols = c("x", "y"), cutoff = 0.1)

# use same mesh for all simulations

#
grid$trend_est = NA
grid$trend_se = NA
for(i in 385:nrow(grid)) {
  # reset seed so each iteration is using same betas
  set.seed(grid$iter[i])
  betas = runif(time_steps, 0, 5)
  sigmaE = exp(log(0.2) + grid$b_trend[i]*1:time_steps)

  s <- sdmTMB_sim(
    x = x, y = y, mesh = mesh, X = X,
    betas = betas, time_steps = time_steps, rho = grid$rho[i],
    phi = grid$obs_sd[i], range = 0.8, sigma_O = 0, sigma_E = sigmaE,
    seed = grid$iter[i], family = tweedie()
  )

  # we could fit all combinations here, e.g. model with/without AR(1) rho, spatial component, trend.
  # but the goal of this is how well we can recover the trend in the presence of those other things
  f <- paste("observed ~ -1 + ", paste(paste0("time",1:time_steps), collapse="+"))

  m <- sdmTMB(
    data = s, formula = formula(f),
    time = "time", spde = mesh,
    ar1_fields = FALSE, include_spatial = FALSE, epsilon_predictor = "time",
    family=tweedie()
  )

  b_indx = grep("\\bb_epsilon\\b",names(m$sd_report$value))[1]
  grid$trend_est[i] = m$sd_report$value[b_indx]
  grid$trend_se[i] = m$sd_report$sd[b_indx]
}
saveRDS(grid, "validation_simulations/tweedie_results.rds")

library(ggplot2)
grid$obs_sd = as.factor(grid$obs_sd)
# shows distribution of trend estimates generally unbiased
grid$phi = grid$obs_sd
# shows distribution of trend estimates generally unbiased
g1 = ggplot(data = dplyr::filter(grid,rho==0), aes(x=as.factor(b_trend), y=trend_est - b_trend)) +
  geom_boxplot(aes(colour = phi), outlier.colour = NA) +
  theme_bw() +
  ylim(-0.1,0.1) +
  xlab("True trend") +
  ylab("Bias, E[trend] - trend")

# shows variance of trend decreases as magnitude of trend inc and obs error decreases
g2 = ggplot(data = dplyr::filter(grid,rho==0), aes(x=as.factor(b_trend), y=log(trend_se))) +
  geom_boxplot(aes(colour = phi), outlier.colour = NA) +
  theme_bw() +
  ylim(c(-4.3,-2)) +
  xlab("True trend") +
  ylab("Ln (trend se)")

pdf("validation_simulations/sim_tweedie.pdf")
gridExtra::grid.arrange(g1, g2)
dev.off()
