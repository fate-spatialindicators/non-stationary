library(ggplot2)
library(dplyr)
library(viridisLite)


df = readRDS("plots/output_figure_1.rds")
my_col = viridis(1, end=0.8)

pdf("plots/Figure_01_motivating.pdf", width = 7, height=6)

ggplot(df, aes(year, mean)) +
  #geom_ribbon(aes(ymin=low,ymax=hi), alpha = 0.3, fill=my_col) +
  geom_point(aes(year,mean),col=my_col,size=2) +
  geom_line(col=my_col) +
  facet_wrap(~common_name, scale="free_y") +
  theme_bw() +
  xlab("Year") +
  ylab(expression(paste("Spatiotemporal ",sigma))) +
  theme(strip.background =element_rect(fill="white"))

dev.off()
