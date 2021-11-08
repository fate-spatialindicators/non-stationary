library(ggplot2)
library(dplyr)
library(viridisLite)


pdf("plots/Figure_02_trends.pdf")

pred = read.csv("output/pred_dens_summary.csv")
# remove duplicates in pred

names = read.csv("data/name_change.csv") %>%
  dplyr::rename(name = old_name)

pred = dplyr::left_join(pred,names) %>%
  dplyr::mutate(trend_sigma = trend/sigma) %>%
  dplyr::filter(loglinear==TRUE, name != c("sole, deepsea")) %>%
  dplyr::arrange(trend_sigma)
pred = pred[which(duplicated(pred$new_name)==FALSE),]
pred$new_name[which(is.na(pred$new_name))] = "Shortspine thornyhead"
#x$name <- factor(x$name, levels = x$name[order(x$val)])

top5_bottom5 = pred[c(1:5, (nrow(pred)-4):nrow(pred)),]
saveRDS(top5_bottom5,"output/top5_bottom5.rds")

pred$new_name <- factor(pred$new_name, levels = pred$new_name)
# change order of factor levels

my_col = viridis(1, end=0.8)
p1 = ggplot(pred, aes(new_name,trend_sigma)) +
  geom_hline(aes(yintercept=0),alpha=0.3, col="red") +
  geom_pointrange(aes(ymin=(trend-2*trend_se)/sigma, ymax=(trend+2*trend_se)/sigma),
                  size=0.8,alpha=0.6,colour=my_col) +
  scale_colour_viridis() +
  xlab("Species") +
  ylab("Trend in spatiotemporal sd (+/- 2SE)") +
  #coord_cartesian(ylim=c(-0.25,0.25)) +
  coord_flip() +
  theme_bw()

p1
dev.off()

######################################
# This makes original without scaled trends

pdf("plots/Figure_S3_trends.pdf")

pred = read.csv("output/pred_dens_summary.csv")
# remove duplicates in pred

names = read.csv("data/name_change.csv") %>%
  dplyr::rename(name = old_name)

pred = dplyr::left_join(pred,names) %>%
  dplyr::filter(loglinear==TRUE, name != c("sole, deepsea")) %>%
  dplyr::arrange(trend)
pred = pred[which(duplicated(pred$new_name)==FALSE),]
pred$new_name[which(is.na(pred$new_name))] = "Shortspine thornyhead"
#x$name <- factor(x$name, levels = x$name[order(x$val)])

top5_bottom5 = pred[c(1:5, (nrow(pred)-4):nrow(pred)),]
saveRDS(top5_bottom5,"output/top5_bottom5.rds")

pred$new_name <- factor(pred$new_name, levels = pred$new_name)
# change order of factor levels

my_col = viridis(1, end=0.8)
p1 = ggplot(pred, aes(new_name,trend)) +
  geom_hline(aes(yintercept=0),alpha=0.3, col="red") +
  geom_pointrange(aes(ymin=(trend-2*trend_se), ymax=(trend+2*trend_se)),
                  size=0.8,alpha=0.6,colour=my_col) +
  scale_colour_viridis() +
  xlab("Species") +
  ylab("Trend in spatiotemporal sd (+/- 2SE)") +
  #coord_cartesian(ylim=c(-0.25,0.25)) +
  coord_flip() +
  theme_bw()

p1
dev.off()
