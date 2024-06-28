rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)

meta <- fread("data/metadata/parsed_patient_metadata.filt.csv")

irep_parsed <- fread("results/irep_out/irep_results.parsed.tsv")
map_df <- fread("results/irep_out/coverage_results.parsed.tsv")
full_df <- fread("data/metadata/irep/irep_pairs.txt") %>%
  select(run_id = V2, species = V3, taxid = V5) %>%
  filter(run_id %in% meta$run_id)


full_df %>%
  left_join(irep_parsed) %>%
  left_join(map_df) %>% 
  left_join(meta %>% select(run_id, hap_vap_cap)) %>%
  group_by(run_id, hap_vap_cap) %>%
  summarise(n_with_reads = n_distinct(species),
            n_deep = sum(mean_depth >= 5, na.rm = T),
            n_rep = sum(bPTR > 0, na.rm = T)) %>%
  filter(n_deep > 0) %>%
  mutate(prop_rep = n_rep / n_deep) %>%
  left_join(meta %>% select(run_id, hosp_los_hours)) %>%
  ggplot(aes(x = log10(as.numeric(hosp_los_hours)), y = prop_rep, color = hap_vap_cap)) +
  geom_point()
  geom_tile()

  full_df %>%
    left_join(irep_parsed) %>%
    left_join(map_df) %>% 
    left_join(meta %>% select(run_id, hap_vap_cap)) %>%
    select(run_id, hap_vap_cap, taxid, species, bPTR) %>%
    ggplot(aes(x = species, y = run_id, fill = bPTR)) +
    geom_tile() +
    facet_grid(rows = vars(hap_vap_cap), space = "free", scale = "free") +
    scale_fill_gradient(na.value = "white") 
  
  full_df %>%
    left_join(irep_parsed) %>%
    left_join(map_df) %>% 
    left_join(meta %>% select(run_id, total_reads = n_reads, hap_vap_cap)) %>%
    group_by(run_id, hap_vap_cap, total_reads) %>%
    summarise(n_with_reads = n_distinct(species),
              n_deep = sum(mean_depth >= 5, na.rm = T),
              n_rep = sum(bPTR > 0, na.rm = T)) %>%
    ggplot(aes(x = total_reads, y = n_rep)) +
    geom_point()

  rep_parsed <- full_df %>%
    left_join(irep_parsed) %>%
    left_join(map_df) %>%
    left_join(meta %>% select(run_id, total_reads = n_reads, hap_vap_cap)) %>%
    filter(hap_vap_cap != "Water control") %>%
    filter(mean_depth >= 1) %>%
    mutate(bPTR = replace_na(bPTR, 0)) %>%
    mutate(is_replicating = ifelse(bPTR == 0, 0, 1)) %>%
    select(run_id, species, is_replicating) %>%
    separate(species, c("genus"), sep = "\\ ")
  
  
  rep_mat <- rep_parsed %>%
    group_by(run_id, genus) %>%
    summarise(is_replicating = sum(is_replicating) > 0) %>%
    pivot_wider(id_cols = run_id, names_from = genus, values_from = is_replicating) %>%
    mutate(across(everything(), ~replace_na(.x, 0))) %>%
    column_to_rownames("run_id")

  genus_filt <- rep_parsed %>%
    group_by(genus) %>%
    summarise(n = n(),
              n_rep = sum(is_replicating)) %>%
    arrange(desc(n)) %>%
    filter(n_rep >= 5) %>%
    filter(n >= 10)
  
  mat <- cor(rep_mat[, genus_filt$genus])

plot_df <- as.data.frame(mat) %>%
  rownames_to_column("genus1") %>%
  pivot_longer(!genus1, names_to = "genus2", values_to = "corr") %>%
  filter(!is.na(corr)) %>%
  # filter(corr > 0.3| corr < -0.3) %>%
  filter(genus1 != genus2)

plot_df %>%
  ggplot(aes(x = genus1, y = genus2, fill = corr)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3))


rep_parsed %>%
  filter(genus == "Paraburkholderia")

rep_parsed %>%
  filter(genus == "Abiotrophia")
rep_mat

cor(rep_mat)
