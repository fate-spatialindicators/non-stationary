library(dplyr)
library(ggplot2)
library(sdmTMB)
#devtools::install_github("seananderson/ggsidekick")
theme_set(ggsidekick::theme_sleek())

species <- read.csv("survey_data/species_list.csv", fileEncoding="UTF-8-BOM")
names(species) <- tolower(names(species))
species <- rename(species, common_name = common.name)
files <- paste0("output/", sub(" ", "_", species$common_name), "_ar1_range15_sigma10.RData")
models_ll <- purrr::map(files, function(f) {
  cat(f, "\n")
  load(f)
  ad_fit_ll
})
models_null <- purrr::map(files, function(f) {
  cat(f, "\n")
  load(f)
  ad_fit
})

get_b_eps <- function(x) {
  est <- as.list(x$sd_report, "Estimate", report = TRUE)$b_epsilon_logit
  se <- as.list(x$sd_report, "Std. Error", report = TRUE)$b_epsilon_logit
  tibble(
    est_b_eps = est,
    lwr_b_eps = est - 1.96 * se,
    upr_b_eps = est + 1.96 * se
  )
}

get_rho <- function(x) {
  lp <- x$tmb_obj$env$last.par.best
  rho <- as.list(x$tmb_obj$report(lp))$rho
  tibble(
    rho = rho
  )
}

get_tweedie_power <- function(x) {
  est <- as.list(x$sd_report, "Estimate")$thetaf
  se <- as.list(x$sd_report, "Std. Error")$thetaf
  tibble(
    est_tweedie_p = plogis(est) + 1,
    lwr_tweedie_p = plogis(est - 1.96 * se) + 1,
    upr_tweedie_p = plogis(est + 1.96 * se) + 1
  )
}

get_tweedie_phi <- function(x) {
  est <- as.list(x$sd_report, "Estimate")$ln_phi
  se <- as.list(x$sd_report, "Std. Error")$ln_phi
  tibble(
    est_phi = exp(est),
    lwr_phi = exp(est - 1.96 * se),
    upr_phi = exp(est + 1.96 * se)
  )
}

b_eps <- purrr::map_dfr(models_ll, get_b_eps)
tweedie_p <- purrr::map_dfr(models_ll, get_tweedie_power)
tweedie_phi <- purrr::map_dfr(models_ll, get_tweedie_phi)
rho <- purrr::map_dfr(models_ll, get_rho)
rho_null <- purrr::map_dfr(models_null, get_rho)
  names(rho_null) <- "rho_null"
sp_names <- tibble(common_name = species$common_name)
out <- bind_cols(sp_names, b_eps, tweedie_p, tweedie_phi, rho, rho_null)

g1 <- out %>%
  filter(est_b_eps > -0.3, est_b_eps < 0.3) %>%
  ggplot(aes(est_b_eps, y = reorder(common_name, est_b_eps), xmin = lwr_b_eps, xmax = upr_b_eps)) +
  geom_pointrange() +
  geom_vline(xintercept = 0, lty = 2) +
  theme(axis.title.y=element_blank())
g1

g2 <- out %>%
  filter(est_b_eps > -0.3, est_b_eps < 0.3) %>%
  ggplot(aes(est_tweedie_p, est_b_eps)) +
  geom_linerange(aes(ymin = lwr_b_eps, ymax = upr_b_eps), alpha = 0.3) +
  geom_linerange(aes(xmin = lwr_tweedie_p, xmax = upr_tweedie_p), alpha = 0.3) +
  geom_point() +
  ggrepel::geom_text_repel(aes(label = common_name), size = 3, alpha = 0.5)
g2

g3 <- out %>%
  filter(est_b_eps > -0.3, est_b_eps < 0.3) %>%
  ggplot(aes(est_phi, est_b_eps)) +
  geom_linerange(aes(ymin = lwr_b_eps, ymax = upr_b_eps), alpha = 0.3) +
  geom_linerange(aes(xmin = lwr_phi, xmax = upr_phi), alpha = 0.3) +
  geom_point() +
  ggrepel::geom_text_repel(aes(label = common_name), size = 3, alpha = 0.3)
g3

g4 <- out %>%
  ggplot(aes(est_tweedie_p, est_phi)) +
  geom_linerange(aes(xmin = lwr_tweedie_p, xmax = upr_tweedie_p), alpha = 0.3) +
  geom_linerange(aes(ymin = lwr_phi, ymax = upr_phi), alpha = 0.3) +
  geom_point() +
  ggrepel::geom_text_repel(aes(label = common_name), size = 3, alpha = 0.3)
g4

g <- cowplot::plot_grid(g1, g4, g2, g3, nrow = 2L)
ggsave("plots/tweedie_vs_epsilon_b.pdf", width = 9.5, height = 9.5)


#plot b_eps and compare rho between null and loglinear models
g5 <- out %>%
  filter(est_b_eps > -0.3, est_b_eps < 0.3) %>%
  ggplot(aes(rho, y = reorder(common_name, est_b_eps), xmin = 0)) +
  geom_point() +
  theme_bw() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank())
g5

g6 <- out %>%
  filter(est_b_eps > -0.3, est_b_eps < 0.3) %>%
  ggplot(aes(rho_null, y = reorder(common_name, est_b_eps), xmin = 0)) +
  geom_point() +
  theme_bw() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank())
g6

g_null <- cowplot::plot_grid(g1, g5, g6, nrow = 1)
ggsave("plots/rho_epsilon_b.pdf", width = 9.5, height = 5)
