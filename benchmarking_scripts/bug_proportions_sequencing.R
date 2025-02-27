rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)

meta_filt <- fread("data/metadata/parsed_patient_metadata.filt.csv") %>%
  filter(high_microbe_count) %>%
  filter(hap_vap_cap %in% c("HAP", "VAP"))

micro_meta <- fread("data/metadata/parsed_microbiology_results.csv")
df_filt <- fread("results/tax_classification_out/abundance_matrices/RA.S.zeroed.decontam.2.csv") %>%
  filter(run_id %in% meta_filt$run_id)

long_df <- df_filt %>%
  pivot_longer(!run_id, names_to = "taxa", values_to = "rel_a") %>%
  filter(rel_a != 0)

# Get dominant microbes
morsels <- foreach(id = unique(long_df$run_id)) %do% {
  temp <- long_df %>%
    filter(run_id == id) %>%
    arrange(desc(rel_a))
  
  temp %>%
    mutate(dominant = c(T, rep(F, nrow(temp) - 1)))
}

parsed <- bind_rows(morsels)
  # filter(!grepl("Candida|Aspergillus|Nakaseomyces", taxa)) 

count_df <- parsed %>%
  group_by(taxa) %>%
  summarise(n = n_distinct(run_id)) %>%
  arrange(desc(n)) %>%
  mutate(prop = n / nrow(meta_filt)) %>%
  head(30)

plot_df <- parsed %>%
  filter(taxa %in% count_df$taxa) %>%
  group_by(taxa) %>%
  summarise(n = n_distinct(run_id)) %>%
  arrange(desc(n)) %>%
  mutate(prop = n / nrow(meta_filt)) %>%
  mutate(is_pathogen = taxa %in% micro_meta$bugs)

plot_df %>%  
  mutate(taxa = factor(taxa, count_df$taxa)) %>%
  ggplot(aes(x = taxa, y = prop, fill = is_pathogen)) +
  geom_bar(stat = "identity", color = "black") +
  theme_classic() +
  scale_fill_manual(values = c("darkseagreen4", "khaki3")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3,
                                   face = "italic")) +
  labs(x = "Bacterial species", y = "Prop. patients", fill = "Detected via PCR/culture")

ggsave("results/benchmarking_out/bug_proportions.sequencing.pdf", 
       dpi = 600, 
       width = 7, 
       height = 3)

plot_df %>% filter(is_pathogen) %>%
  arrange(desc(prop))
