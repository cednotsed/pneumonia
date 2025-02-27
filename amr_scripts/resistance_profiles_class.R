rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(ggpubr)
require(viridis)
require(Hmisc)

meta <- fread("data/metadata/parsed_patient_metadata.filt.csv") %>%
  filter(high_microbe_count)

df <- fread("results/amr_out/amr_matrices/resfinder.class.decontam.csv") %>%
  pivot_longer(!run_id, names_to = "class", values_to = "presence") %>%
  inner_join(meta %>% select(run_id, hap_vap_cap, hap_vap2)) %>%
  filter(hap_vap_cap %in% c("HAP", "VAP"))

count_df <- df %>%
  group_by(hap_vap2) %>%
  summarise(n_total = n_distinct(run_id))

plot_df <- df %>%
  group_by(hap_vap2, class) %>%
  summarise(n = sum(presence == 1)) %>%
  arrange(desc(n)) %>%
  left_join(count_df) %>%
  mutate(prop = n / n_total)

order_df <- plot_df %>%
  filter(hap_vap2 == "VAP") %>%
  arrange(desc(prop))

plot_df %>%  
  mutate(class = factor(class, unique(order_df$class))) %>%
  ggplot(aes(x = class, y = hap_vap2, fill = prop)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1)) +
  labs(x = "Predicted resistance", y = "HAP subtype", fill = "Prop. patients") 

ggsave("results/amr_out/amr_profiles.pdf", width = 9, height = 3, dpi = 600)