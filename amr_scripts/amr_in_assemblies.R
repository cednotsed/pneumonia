rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(ggpubr)
require(viridis)
require(Hmisc)

meta <- fread("data/metadata/parsed_patient_metadata.filt.csv") %>%
  filter(high_microbe_count) %>%
  filter(hap_vap_cap %in% c("HAP", "VAP"))

genome_meta <- fread("results/assembly_out/assembly_merged_metadata.csv") %>%
  filter(gtdbtk_warnings == "N/A") %>%
  filter(run_id %in% meta$run_id) %>%
  filter(checkm2_contamination <= 5) %>%
  filter(checkm2_completeness >= 90) %>%
  filter(species != "")

genome_meta %>%
  group_by(species) %>%
  summarise(n = n()) %>%
  arrange(desc(n))

amr_df <- fread("results/amr_out/amr_matrices/resfinder.assemblies.gene_content.csv") %>%
  filter(genome_name %in% genome_meta$genome_name)

# Median classes
sum_df <- amr_df %>%
  group_by(genome_name) %>%
  summarise(n = n_distinct(class)) %>%
  ungroup() %>%
  arrange(desc(n))

parsed <- genome_meta %>%
  select(genome_name) %>%
  left_join(sum_df) %>%
  mutate(n = replace_na(n, 0))

parsed %>%
  summarise(max = max(n),
            median = median(n),
            prop_multi = sum(n > 1) / nrow(parsed))

parsed %>%
  mutate(multi = ifelse(n > 1, "Yes", "No")) %>%
  ggplot(aes(x = n, fill = multi)) +
  geom_histogram(color = "black") +
  scale_fill_manual(values = c("grey", "indianred")) +
  theme_bw() +
  labs(x = "No. predicted resistance classes", y = "No. patients", fill = "Is multi-drug resistant?")


ggsave("results/amr_out/assembly_AMR_classes_histogram.pdf", dpi = 600, width = 5, height = 4)

amr_sum_df <- merged %>%
  group_by(genome_name, species, antimicrobial) %>%
  summarise(resistant = ifelse(sum(wgs_predicted_phenotype == "Resistant") > 0, 1, 0)) %>%
  group_by(genome_name, species) %>%
  summarise(n = sum(resistant)) %>%
  ungroup() %>%
  arrange(desc(n))

amr_sum_df %>%
  ggplot(aes(x = species, y = n)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3))
summarise(max = max(n),
            median = median(n))
