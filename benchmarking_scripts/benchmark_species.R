rm(list = ls())
setwd("C:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)

micro_meta <- fread("data/metadata/parsed_microbiology_results.bacterial_sp_only.csv") %>%
  filter(species_resolved)

meta_filt <- fread("data/metadata/parsed_patient_metadata.filt.csv") %>%
  filter(hap_vap_cap %in% c("HAP", "VAP")) %>%
  filter(high_microbe_count)

df_filt <- fread("results/tax_classification_out/abundance_matrices/RA.S.relaxed.csv") %>%
  filter(run_id %in% meta_filt$run_id) %>%
  filter(run_id %in% micro_meta$run_id)

long_df <- df_filt %>%
  pivot_longer(!run_id, names_to = "taxa", values_to = "rel_a") %>%
  left_join(meta_filt) %>%
  filter(rel_a != 0)

morsels <- foreach(id = df_filt$run_id) %do% {
  # id = "27_4"
  micro_temp <- micro_meta %>%
    filter(run_id == id)
  
  temp <- long_df %>%
    filter(run_id == id)
  
  total_df <- micro_temp %>%
    group_by(method) %>%
    summarise(n_identified = sum(bugs != "Negative"))
  # bind_rows(tibble(method = "sequencing", n_identified = n_distinct(temp$genus)))
  
  itx_df <- micro_temp %>%
    group_by(method) %>%
    summarise(n_intersect = sum(bugs %in% temp$taxa)) %>%
    left_join(total_df) %>%
    mutate(run_id = id)
  
  return(itx_df)
}

bind_rows(morsels) %>% 
  # filter(method == "culture") %>%
  # filter(n_identified != 0, n_intersect < n_identified)
  group_by(method) %>%
  summarise(sum_sequenced = sum(n_intersect),
            sum_total = sum(n_identified)) %>%
  mutate(sum_sequenced / sum_total)




