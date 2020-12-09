library(sdmTMB)
library(sp)
library(ggplot2)
library(ggsidekick)
species = read.csv("survey_data/species_list.csv")
bio = readRDS("survey_data/wcbts_bio_2019-08-01.rds")
haul = readRDS("survey_data/wcbts_haul_2019-08-01.rds")
catch = readRDS("survey_data/wcbts_catch_2019-08-01.rds")
names(catch) = tolower(names(catch))

catch$trawl_id = as.character(catch$trawl_id)
haul$trawl_id = as.character(haul$trawl_id)
dat = dplyr::left_join(catch, haul)
# filter species from list
dat = dplyr::filter(dat, scientific_name %in% species$Scientific.Name)
# convert to UTM
coordinates(dat) <- c("longitude_dd", "latitude_dd")
proj4string(dat) <- CRS("+proj=longlat +datum=WGS84")
newproj = paste("+proj=utm +zone=10 ellps=WGS84")
dat <- spTransform(dat, CRS(newproj))
dat = as.data.frame(dat)
dat$lon = dat$longitude_dd/1000
dat$lat = dat$latitude_dd/1000
#dat$year = as.numeric(substr(dat$date_yyyymmdd,1,4))

dat = dplyr::filter(dat, scientific_name %in% c("Eopsetta jordani",
  "Sebastolobus altivelis","Sebastes crameri","Sebastes jordani"))

models = list()
for(i in 1:length(unique(dat$scientific_name))) {
  
  sub = dplyr::filter(dat, scientific_name==unique(dat$scientific_name)[i]) %>% 
    dplyr::filter(cpue_kg_km2>0)
  spde <- make_mesh(sub, c("lon", "lat"), cutoff = 10) 
  sub$depth_scaled = as.numeric(scale(sub$depth_m))
  sub$depth_scaled2 = sub$depth_scaled^2

  m = sdmTMB(cpue_kg_km2 ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
    data = sub, time = "year", spde = spde, family = Gamma(link = "log"))
  models[[i]] = m
  
  # calculate residuals 
  sub$resid = sub$cpue_kg_km2 - predict(m)$est
  # bootstrap residuals by year to calculate sd broom::bootstrap() deprecated
  sp_df = data.frame("species"=sub$scientific_name[1], "year"=unique(dat$year), "low"=NA, "mean"=NA, "hi"=NA)
  for(y in 1:length(unique(dat$year))) {
    sub_yr = dplyr::filter(sub, year==unique(dat$year)[y]) %>% dplyr::select(resid)
    m = matrix(sample(sub_yr$resid, size = length(sub_yr$resid)*1000, replace=TRUE), 1000, length(sub_yr$resid))
    sds = apply(m,1,sd)
    sp_df$low[y] = quantile(log(sds),0.025)
    sp_df$mean[y] = mean(log(sds))
    sp_df$hi[y] = quantile(log(sds),0.975)
  }
  if(i==1) {
    df = sp_df
  } else {
    df = rbind(df, sp_df)
  }
}

# bring in common names
df$common_name = ""
df$common_name[which(df$species=="Eopsetta jordani")] = "Petrale sole"
df$common_name[which(df$species=="Sebastolobus altivelis")] = "Longspine thornyhead"
df$common_name[which(df$species=="Sebastes crameri")] = "Darkblotched rockfish"
df$common_name[which(df$species=="Sebastes jordani")] = "Shortbelly rockfish"


g1 = df %>%
ggplot(aes(year, mean)) + 
  geom_ribbon(aes(ymin=low,ymax=hi),alpha=0.3,fill="darkblue") + 
  geom_point(aes(year,mean),col="darkblue",size=2) + 
  geom_line(col="darkblue") + 
  facet_wrap(~common_name,scale="free_y") + 
  theme_sleek() + 
  xlab("Year") + ylab(expression(paste("Ln(",sigma,") of residuals")))

pdf("plots/figure_1.pdf") 
g1
dev.off()

jpeg("plots/figure_1.jpeg") 
g1
dev.off()
