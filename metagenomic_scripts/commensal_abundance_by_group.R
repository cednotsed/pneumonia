rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(randomcoloR)

meta_filt <- fread("data/metadata/parsed_patient_metadata.filt.csv") %>%
  filter(high_microbe_count) %>%
  mutate(sample_type = ifelse(hap_vap_cap == "healthy", "SPU", sample_type)) %>%
  filter(sample_type == "SPU"|hap_vap_cap == "Water control")

df_filt <- fread("results/tax_classification_out/abundance_matrices/RA.G.zeroed.csv") %>%
  filter(run_id %in% meta_filt$run_id)

long_df <- df_filt %>%
  pivot_longer(!run_id, names_to = "taxa", values_to = "rel_a") %>%
  left_join(meta_filt)

plot_df <- long_df %>% 
  filter(taxa %in% c("Prevotella", "Streptococcus", "Veillonella", 
                     "Fusobacterium", "Haemophilus", "Veillonella", 
                     "Rothia")) %>%
  mutate(hap_vap_cap = ifelse(hap_vap_cap == "healthy", 
                              "Non-pneumonia control", 
                              hap_vap_cap))

count_df <- plot_df %>%
  group_by(hap_vap_cap) %>%
  summarise(n = n_distinct(run_id))

order_df <- plot_df %>%
  group_by(run_id, hap_vap_cap) %>%
  summarise(sum_rel_a = sum(rel_a)) %>%
  arrange(desc(sum_rel_a))
 
order_df %>%
  mutate(hap_vap_cap = factor(hap_vap_cap, c("HAP", "VAP", "CAP", 
                                             "Non-pneumonia control", "Water control"))) %>%
  ggplot(aes(x = hap_vap_cap, y = sum_rel_a, fill = hap_vap_cap)) +
  geom_boxplot() +
  geom_text(aes(y = 1.25, label = str_glue("n={n}")),
            data = count_df) +
  theme_bw() +
  scale_fill_manual(values = list(HAP = "blue", VAP = "goldenrod", 
                                  CAP = "red", `Non-pneumonia control` = "olivedrab", 
                                  `Water control` = "grey44")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Sample type", y = "Commensal relative abundance") +
  theme(legend.position = "none")
 
long_df %>%
  filter(hap_vap_cap == "healthy") %>%
  filter(taxa == "Aspergillus") %>%
  filter(rel_a != 0)
# plot_df %>%
#   mutate(run_id = factor(run_id, order_df$run_id)) %>%
#   ggplot(aes(x = run_id, y = rel_a, fill = taxa)) +
#   geom_bar(stat = "identity") +
#   facet_grid(cols = vars(hap_vap_cap), space = "free", scale = "free")

