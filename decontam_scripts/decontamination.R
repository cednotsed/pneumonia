rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(vegan)
require(ape)

meta <- fread("data/metadata/parsed_patient_metadata.filt.csv") %>%
  filter(high_microbe_count) %>%
  filter(run != 1) %>%
  full_join(tibble(run_id = "A_7", run = c("A", "B", "C", "D"),
                   hap_vap_cap = "Water control"))

# Parse sequencing results
patient_meta <- meta %>% 
  filter(hap_vap_cap %in% c("HAP", "VAP", "CAP", "healthy"))

control_meta <- meta %>%
  filter(hap_vap_cap == "Water control")

rank <- "G"
threshold <- 2
RA <- fread(str_glue("results/tax_classification_out/abundance_matrices/RA.{rank}.zeroed.csv"))
read <- fread(str_glue("results/tax_classification_out/abundance_matrices/read_counts.{rank}.zeroed.csv"))

RA_filt <- RA %>%
  filter(run_id %in% patient_meta$run_id)

read_filt <- read %>%
  filter(run_id %in% patient_meta$run_id)

# Use run A for all healthy runs
control_RA <- RA %>%
  filter(run_id %in% control_meta$run_id)

control_read <- read %>%
  filter(run_id %in% control_meta$run_id)

long_RA <- RA_filt %>%
  pivot_longer(!run_id, names_to = "taxa", values_to = "rel_a")

long_read <- read_filt %>%
  pivot_longer(!run_id, names_to = "taxa", values_to = "read_count")

control_long_RA <- control_RA %>%
  pivot_longer(!run_id, names_to = "taxa", values_to = "control_rel_a") %>%
  left_join(control_meta %>% distinct(run_id, run)) %>% 
  select(taxa, run, control_rel_a) %>%
  group_by(taxa, run) %>%
  summarise(control_rel_a = max(control_rel_a))

decontam_df <- long_RA %>%
  left_join(meta %>% select(run_id, run)) %>%
  filter(rel_a != 0) %>%
  left_join(control_long_RA) %>%
  left_join(long_read) %>%
  select(run_id, taxa, rel_a, read_count, control_rel_a) %>% 
  filter(rel_a > threshold * control_rel_a)

View(decontam_df)
decontam_read <- decontam_df %>%
  select(-rel_a, -control_rel_a) %>%
  pivot_wider(id_cols = run_id, names_from = taxa, values_from = read_count) %>%
  mutate(across(everything(), ~replace_na(., 0)))

# Rescale relative abundance
otu_to_RA <- function(df) {
  mat <- as.matrix(df)
  RA_df <- as.data.frame(mat / rowSums(mat))
  colnames(RA_df) <- colnames(df)
  
  return(RA_df)
}

decontam_RA <- otu_to_RA(decontam_read %>% column_to_rownames("run_id")) %>%
  rownames_to_column("run_id")

decontam_RA %>%
  fwrite(str_glue("results/tax_classification_out/abundance_matrices/RA.{rank}.zeroed.decontam.{threshold}.csv"))

decontam_read %>%
  fwrite(str_glue("results/tax_classification_out/abundance_matrices/read_counts.{rank}.zeroed.decontam.{threshold}.csv"))

