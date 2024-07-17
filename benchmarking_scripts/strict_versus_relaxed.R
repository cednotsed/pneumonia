rm(list = ls())
setwd("C:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)

meta_filt <- fread("data/metadata/parsed_patient_metadata.filt.csv") %>%
  filter(hap_vap_cap %in% c("HAP", "VAP", "CAP")) %>%
  filter(high_microbe_count)

relaxed_df <- fread("results/tax_classification_out/abundance_matrices/RA.G.relaxed.csv") %>%
  filter(run_id %in% meta_filt$run_id) %>%
  pivot_longer(!run_id, names_to = "taxa", values_to = "rel_a") %>%
  left_join(meta_filt) %>%
  filter(rel_a != 0)

strict_df <- fread("results/tax_classification_out/abundance_matrices/RA.G.zeroed.decontam.2.csv") %>%
  filter(run_id %in% meta_filt$run_id) %>%
  pivot_longer(!run_id, names_to = "taxa", values_to = "rel_a") %>%
  left_join(meta_filt) %>%
  filter(rel_a != 0)

relaxed_df %>%
  group_by(run_id) %>%
  summarise(n = n_distinct(taxa)) %>%
  ungroup() %>%
  summarise(median = median(n),
            lower_q = quantile(n, 0.25),
            higher_q = quantile(n, 0.75))

strict_df %>%
  group_by(run_id) %>%
  summarise(n = n_distinct(taxa)) %>%
  ungroup() %>%
  summarise(median = median(n),
            lower_q = quantile(n, 0.25),
            higher_q = quantile(n, 0.75))

