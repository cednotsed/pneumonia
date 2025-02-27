rm(list = ls())
setwd("C:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(ggpubr)
require(viridis)

meta <- fread("data/metadata/parsed_patient_metadata.filt.csv") %>%
  filter(hap_vap_cap %in% c("HAP", "VAP"))

meta %>%
  filter(hap_vap2 == "NV-HAP", sample_type == "ND-BAL")

plot_df <- meta %>%
  group_by(hap_vap2, sample_type) %>%
  summarise(n = n())
sample_counts <- meta %>%
  group_by(hap_vap2) %>%
  summarise(n_samples = n())

plot_df %>%  
  group_by(hap_vap2) %>%
  mutate(prop = n / sum(n)) %>%
  ggplot(aes(x = hap_vap2, y = prop)) +
    geom_bar(aes(fill = sample_type), stat = "identity", color = "black") +
    theme_bw() +
  labs(x = "HAP subtype", y = "% samples") +
  geom_text(aes(x = hap_vap2, y = 1.05, label = str_glue("n={n_samples}")),
            data = sample_counts)

ggsave("results/summary_out/sample_type_distribution.pdf", width = 5, height = 4, dpi = 600)

mat <- meta %>%
  mutate(sample_type = ifelse(sample_type == "SPU", "SPU", "non-SPU")) %>%
  group_by(ventilation, sample_type) %>%
  summarise(n = n()) %>%
  pivot_wider(id_cols = ventilation, names_from = sample_type, values_from = n) %>%
  column_to_rownames("ventilation")

test <- fisher.test(mat)

mat["SPU"] / apply(mat, 1, sum)


