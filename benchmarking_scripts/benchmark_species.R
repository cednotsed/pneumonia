rm(list = ls())
setwd("C:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)

to_assess <- fread("results/benchmarking_out/to_assess_ids.txt")

micro_meta <- fread("data/metadata/parsed_microbiology_results.bacterial_sp_only.csv") %>%
  filter(species_resolved) %>%
  filter(run_id %in% to_assess$run_id)

n_distinct(micro_meta$run_id)

meta_filt <- fread("data/metadata/parsed_patient_metadata.filt.csv") %>%
  filter(hap_vap_cap %in% c("HAP", "VAP")) %>%
  filter(high_microbe_count) %>%
  filter(run_id %in% to_assess$run_id)

n_distinct(meta_filt$run_id)

df_filt <- fread("results/tax_classification_out/abundance_matrices/RA.S.zeroed.2.csv") %>%
  filter(run_id %in% meta_filt$run_id) %>%
  filter(run_id %in% micro_meta$run_id)

long_df <- df_filt %>%
  pivot_longer(!run_id, names_to = "taxa", values_to = "rel_a") %>%
  left_join(meta_filt) %>%
  filter(rel_a != 0)

morsels <- foreach(id = df_filt$run_id) %do% {
  # id = "1_3"
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

plot_df <- bind_rows(morsels) %>% 
  mutate(prop_identified = n_intersect / n_identified)

plot_df %>%
  group_by(method) %>%
  summarise(n_poor = sum(n_intersect == n_identified) / n_distinct(run_id),
            n_patients = n_distinct(run_id))

plot_df
ggplot(aes(x = method, y = prop_identified)) +
  geom_boxplot()
  group_by(method) %>%
  summarise(median = median(prop_identified))

plot_df %>%
  group_by(method) %>%
  summarise(sum_sequenced = sum(n_intersect),
            sum_total = sum(n_identified)) %>%
  mutate(sum_sequenced / sum_total)




