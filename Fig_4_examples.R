library(sdmTMB)
library(dplyr)
library(ggplot2)
library(viridis)


# Species of interest
species = read.csv("survey_data/species_list.csv", fileEncoding="UTF-8-BOM")
names(species) = tolower(names(species))
species = dplyr::rename(species,
                        common_name = common.name,
                        scientific_name = scientific.name)

# all predictions
d = readRDS("output/predictions_all.rds")

# filter out a few species with trends != 0
sub = dplyr::filter(d, common_name %in% c("sole, rex"), lat > 4500)

# version 1: raw epsilon_st values
g0 = sub %>%
  dplyr::filter(year%in%c(2003,2018)) %>%
  ggplot(aes(lon,lat,fill=epsilon_st)) +
  geom_tile() +
  scale_fill_gradient2() +
  facet_grid(model~year) +
  theme_bw() +
  theme(strip.background =element_rect(fill="white"),
        strip.text.x = element_text(size = 8)) +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle("Raw epsilon_st")

# option 2: show percentiles
summary_stats = dplyr::group_by(sub, common_name) %>%
  dplyr::summarize(mu = mean(epsilon_st),
                   sd = sd(epsilon_st))
sub = dplyr::left_join(sub, summary_stats)
sub$percentile = pnorm(sub$epsilon_st, sub$mu, sub$sd)

g1 = sub %>%
  dplyr::filter(year%in%c(2003,2018)) %>%
  ggplot(aes(lon,lat,fill=percentile)) +
  geom_tile() +
  scale_fill_gradient2() +
  facet_grid(model~year) +
  theme_bw() +
  theme(strip.background =element_rect(fill="white"),
        strip.text.x = element_text(size = 8)) +
  xlab("Longitude") + ylab("Latitude") + ggtitle("Percentile")

# Same plot as g1 but log percentile
g2 = sub %>%
  dplyr::filter(year%in%c(2003,2018)) %>%
  ggplot(aes(lon,lat,fill=log(percentile))) +
  geom_tile() +
  scale_fill_gradient2() +
  facet_grid(model~year) +
  theme_bw() +
  theme(strip.background =element_rect(fill="white"),
        strip.text.x = element_text(size = 8)) +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle("Log percentile")

summary_stats2 = sub %>%
  dplyr::filter(year==2003) %>%
  group_by(common_name, model) %>%
  dplyr::summarize(sd_0 = sd(epsilon_st)) #%>%
  #dplyr::select(common_name,sd_0)
sub = dplyr::left_join(sub, summary_stats2)
sub$scaled_epsilon_st = sub$epsilon_st/sub$sd_0

g3 = sub %>%
  dplyr::filter(year%in%c(2003,2018)) %>%
  ggplot(aes(lon,lat,fill=scaled_epsilon_st)) +
  geom_tile() +
  scale_fill_gradient2() +
  facet_grid(model~year) +
  theme_bw() +
  theme(strip.background =element_rect(fill="white"),
        strip.text.x = element_text(size = 8)) +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle("Standardized to epsilon_2003")

sub = dplyr::group_by(sub,cell) %>%
  dplyr::mutate(demeaned = epsilon_st - mean(epsilon_st))

g4 = sub %>%
  dplyr::filter(year%in%c(2003,2018)) %>%
  ggplot(aes(lon,lat,fill=demeaned)) +
  geom_tile() +
  scale_fill_gradient2() +
  facet_grid(model~year) +
  theme_bw() +
  theme(strip.background =element_rect(fill="white"),
        strip.text.x = element_text(size = 8)) +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle("Demeaned across years/models")

pdf("Fig_4_examples_lingcod.pdf")
g0
g1
g2
g3
g4
dev.off()
