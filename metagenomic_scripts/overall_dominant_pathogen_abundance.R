rm(list = ls())
setwd("C:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(ggpubr)
require(viridis)

meta <- fread("data/metadata/parsed_patient_metadata.filt.csv")

meta_filt <- meta %>%
  filter(hap_vap_cap %in% c("HAP", "VAP")) %>%
  filter(high_microbe_count) %>%
  mutate(ventilation = ifelse(ventilation, "VAP", "HAP"))

df_filt <- fread("results/tax_classification_out/abundance_matrices/RA.G.zeroed.decontam.2.csv") %>%
  filter(run_id %in% meta_filt$run_id)

long_df <- df_filt %>%
  pivot_longer(!run_id, names_to = "taxa", values_to = "rel_a") %>%
  filter(rel_a != 0) %>%
  left_join(meta_filt %>% select(run_id, ventilation)) 

plot_df <- long_df %>%
  group_by(run_id) %>%
  summarise(max_a = max(rel_a))
  
plot_df %>% summarise(median = median(max_a),
                      lower = quantile(max_a, 0.25),
                      upper = quantile(max_a, 0.75))
plot_df %>%
ggplot(aes(x = 1, y = max_a)) +
  geom_violin() +
  geom_boxplot(width = 0.3)

               