rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)

read_counts <- fread("results/qc_out/sample_read_counts.csv")

meta <- fread("data/metadata/parsed_patient_metadata.csv") %>%
  filter(!grepl("None", run_id)) %>%
  left_join(read_counts) %>%
  mutate(n_reads = ifelse(is.na(n_reads), 0, n_reads))

dups <- meta %>%
  filter(duplicated(sample_id) & 
           !grepl("water", sample_id, ignore.case = T))

dup_df <- meta %>%
  filter(sample_id %in% dups$sample_id)

non_dup_df <- meta %>%
  filter(!(sample_id %in% dups$sample_id))

dup_filt <- foreach(sample_name = unique(dup_df$sample_id), .combine = "bind_rows") %do% {
  dup_df %>%
    filter(sample_id == sample_name) %>%
    arrange(desc(n_reads)) %>%
    head(1)
}

merged <- bind_rows(non_dup_df, dup_filt)

# Deduplicate water controls
non_controls <- merged %>%
  filter(hap_vap_cap != "Water control")

water_controls <- merged %>%
  filter(hap_vap_cap == "Water control")

water_control_dedup <- foreach(run_name = unique(water_controls$run)) %do% {
  temp <- water_controls %>%
    filter(run == run_name) %>%
    arrange(desc(n_reads)) %>%
    head(1)
}

final <- bind_rows(non_controls, water_control_dedup)

# Get total and microbial read counts
dat <- fread(str_glue("results/tax_classification_out/abundance_matrices/abundance_matrix.G.tsv")) %>%
  select(-any_of(c("Homo sapiens", "Homo")), -unclassified) %>%
  as_tibble() %>%
  column_to_rownames("run_id") 

microbe_df <- tibble(run_id = rownames(dat),
                     microbial_reads = rowSums(dat))

final %>%
  left_join(microbe_df) %>% 
  mutate(microbial_reads = ifelse(is.na(microbial_reads), 0, microbial_reads)) %>%
  mutate(high_microbe_count = microbial_reads >= 100) %>%
  fwrite("data/metadata/parsed_patient_metadata.filt.csv")

