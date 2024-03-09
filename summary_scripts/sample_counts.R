setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(doParallel)
require(ggforce)
require(ggpubr)
require(scales)
require(ggsci)

meta <- fread("data/metadata/parsed_patient_metadata.csv")

plot_df <- meta %>%
  group_by(hap_vap_cap) %>%
  summarise(cnt = n(), total = nrow(meta))

rpie <- 0.7
rlabel <-  0.6 * rpie

plot_df %>%
  mutate(end_angle = 2 * pi * cumsum(cnt) / total,      # ending angle for each pie slice
         start_angle = lag(end_angle, default = 0),   # starting angle for each pie slice
         mid_angle = 0.5 * (start_angle + end_angle)) %>%  # middle of each pie slice, for the text label
  ggplot() +
  geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = rpie,
                   start = start_angle, end = end_angle, fill = hap_vap_cap)) +
  geom_text(aes(x = rlabel * sin(mid_angle), y = rlabel * cos(mid_angle), label = cnt),
            hjust = 0.5, vjust = 1, size = 5) +
  coord_fixed() +
  scale_x_continuous(limits = c(-1, 1), name = "", breaks = NULL, labels = NULL) +
  scale_y_continuous(limits = c(-1, 1), name = "", breaks = NULL, labels = NULL) +
  scale_fill_manual(values = c("firebrick2", "mediumorchid3", "dodgerblue3", "grey90")) +
  theme(panel.background = element_rect(fill = "white"),
        plot.margin=grid::unit(c(0,0,0,0),"cm"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks=element_blank(),
        axis.ticks.length = unit(0, "pt")) +
  labs(x = NULL, y = NULL, fill = "Sample type")

ggsave("results/metagenomic_out/sample_pie.png",
       width = 5, height = 5)  

