rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(vegan)
require(ape)

meta_filt <- fread("data/metadata/parsed_patient_metadata.filt.csv") %>%
  filter(high_microbe_count) %>%
  filter(hap_vap_cap %in% c("HAP", "VAP")) 

RA_filt <- fread("results/tax_classification_out/abundance_matrices/RA.G.zeroed.decontam.2.csv") %>%
  inner_join(meta_filt %>% select(run_id, hap_vap_cap)) %>%
  select(-hap_vap_cap) %>%
  column_to_rownames("run_id")

dim(RA_filt)

# Bray-Curtis PCoA
bc <- vegdist(as.matrix(RA_filt), method = "bray")

PCOA <- pcoa(bc)

eigenvals <- round(PCOA$values$Relative_eig * 100, 1)

plot_df <- as.data.frame(PCOA$vectors) %>%
  rownames_to_column("run_id") %>%
  left_join(meta_filt)

min(plot_df$Axis.1)
max(plot_df$Axis.1)
min(plot_df$Axis.2)
max(plot_df$Axis.2)

p1 <- plot_df %>%
  ggplot(aes(Axis.1, Axis.2, color = ventilation, shape = ventilation)) +
  geom_point(alpha = 0.8, size = 3) +
  # scale_color_manual(values = c("goldenrod", "blue")) +
  theme_classic() + 
  ylim(-0.6, 0.7) +
  xlim(-0.7, 0.7) +
  labs(x = str_glue("PCo1 ({eigenvals[1]}%)"),
       y = str_glue("PCo2 ({eigenvals[2]}%)")) +
  theme(legend.position = "none",
        axis.title = element_blank())

a1 <- plot_df %>%
  ggplot(aes(y = Axis.1, x = ventilation, fill = ventilation)) +
  geom_boxplot() +
  # scale_fill_manual(values = c("goldenrod", "blue")) +
  theme_bw() +
  geom_pwc(label.size = 3, step.increase = 0.05) +
  coord_flip() +
  ylim(-0.7, 0.7) +
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_blank()) +
  labs(y = str_glue("PCo1 ({eigenvals[1]}%)"))

a2 <- plot_df %>%
  ggplot(aes(x = ventilation, y = Axis.2, fill = ventilation)) +
  geom_boxplot() +
  # scale_fill_manual(values = c("goldenrod", "blue")) +
  theme_bw() +
  geom_pwc(label.size = 3, step.increase = 0.05) +
  ylim(-0.6, 0.7) +
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank()) +
  labs(y = str_glue("PCo2 ({eigenvals[2]}%)"))

ggarrange(a2, p1, ggplot(), a1, nrow = 2, ncol = 2,
          heights = c(3, 1), widths = c(1, 4),
          align = "hv")
