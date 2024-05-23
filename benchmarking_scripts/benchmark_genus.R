rm(list = ls())
setwd("C:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)

micro_meta <- fread("data/metadata/parsed_microbiology_results.csv") %>%
  filter(species_resolved)

meta_filt <- fread("data/metadata/parsed_patient_metadata.filt.csv")

df_filt <- fread("results/metagenomic_out/RA.S.zeroed.csv") %>%
  filter(run_id %in% meta_filt$run_id) %>%
  filter(run_id %in% micro_meta$run_id)

long_df <- df_filt %>%
  pivot_longer(!run_id, names_to = "genus", values_to = "rel_a") %>%
  left_join(meta_filt) 

long_df %>%
  filter(grepl("Candida|Aspergillus", genus)) %>% 
  filter(rel_a != 0) %>% View()
morsels <- foreach(id = df_filt$run_id) %do% {
  # id = "7_5"
  micro_temp <- micro_meta %>%
    filter(run_id == id)
  
  temp <- long_df %>%
    filter(run_id == id)
  
  total_df <- micro_temp %>%
    group_by(method) %>%
    summarise(n_identified = sum(bug_genus != "Negative"))
    # bind_rows(tibble(method = "sequencing", n_identified = n_distinct(temp$genus)))
  
  itx_df <- micro_temp %>%
    group_by(method) %>%
    summarise(n_intersect = sum(bug_genus %in% temp$genus)) %>%
    left_join(total_df) %>%
    mutate(run_id = id)
  
  return(itx_df)
}

bind_rows(morsels) %>% 
  # filter(method == "curetis") %>%
  # filter(n_identified != 0, n_intersect < n_identified)
  group_by(method) %>%
  summarise(sum_sequenced = sum(n_intersect),
            sum_total = sum(n_identified)) %>%
  mutate(sum_sequenced / sum_total)

bind_rows(morsels) %>% 
  filter(method == "biofire") %>%
  filter(n_identified != 0, n_intersect < n_identified)

max(genus_reads)

meta <- fread("results/qc_out/high_microbe_samples.G.csv")

meta %>% View()
hist(log10(meta$microbial_reads))
micro_meta %>% 
  distinct(bugs) %>% View()

micro_meta %>%
  distinct(bugs)

micro
