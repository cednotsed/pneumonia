rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(vegan)
require(hillR)
require(ggpubr)
require(randomcoloR)

meta <- fread("data/metadata/parsed_patient_metadata.filt.csv") %>%
  filter(high_microbe_count)

file_dir <- "results/amr_out/resfinder_out/all_resfinder_reports/"
file_list <- list.files(file_dir, full.names = T)

res_df <- foreach(file = file_list, .combine = "bind_rows") %do% {
  id <- gsub(file_dir, "", file)
  
  fread(file) %>%
    mutate(run_name = id)
}

parsed <- res_df %>%
  rename_all(~tolower(gsub(" |-", "_", .x))) %>%
  dplyr::rename(antimicrobial = `#_antimicrobial`) %>%
  mutate(run_id = gsub("INHALE_FRESH_|.no_human|Library|_raw", "", run_name)) %>%
  mutate(run_id = gsub(".pheno_table.txt", "", run_id)) %>%
  mutate(run_id = gsub("barcode", "", run_id)) %>%
  mutate(run_id = gsub("-0", "-", run_id)) %>%
  mutate(run_id = gsub("-", "_", run_id)) %>%
  mutate(run_id = gsub("12a", "12", run_id)) %>%
  inner_join(meta %>% select(run_id, hap_vap_cap, run)) %>%
  select(-run_name) %>%
  arrange(desc(antimicrobial))


parsed %>%
  fwrite("results/amr_out/resfinder_out/resfinder.parsed.csv")

# Decontaminate
parsed_filt <- parsed %>%
  filter(run != 1) %>%
  filter(hap_vap_cap != "healthy")

morsels <- foreach(run_name = unique(parsed_filt$run)) %do% {
  controls <- parsed_filt %>%
    filter(hap_vap_cap == "Water control",
           run == run_name) %>%
    filter(wgs_predicted_phenotype == "Resistant")
  
  patients <- parsed_filt %>%
    filter(hap_vap_cap != "Water control",
           run == run_name)
  
  decontam <- patients %>%
    filter(!(antimicrobial %in% controls$antimicrobial))
  
  return(decontam)
}

decontam_df <- bind_rows(morsels)

decontam_df %>%
  fwrite("results/amr_out/resfinder_out/resfinder.parsed.decontam.csv")

# Antibiotic matrix
mat <- decontam_df %>%
  select(run_id, antimicrobial, wgs_predicted_phenotype) %>%
  mutate(wgs_predicted_phenotype = ifelse(wgs_predicted_phenotype == "Resistant", 1, 0)) %>%
  mutate(antimicrobial = gsub("\\+|\\-", "_", antimicrobial)) %>%
  mutate(antimicrobial = gsub("\\ ", "_", antimicrobial)) %>%
  pivot_wider(id_cols = run_id, names_from = antimicrobial, values_from = wgs_predicted_phenotype)

mat[is.na(mat)] <- 0

mat %>%
  fwrite("results/amr_out/amr_matrices/resfinder.amr.decontam.csv")

# Resistance class matrix
mat <- decontam_df %>%
  group_by(run_id, class) %>%
  summarise(wgs_predicted_phenotype = ifelse(sum(wgs_predicted_phenotype == "Resistant", na.rm = T) > 0,
                                             1, 0)) %>%
  pivot_wider(id_cols = run_id, names_from = class, values_from = wgs_predicted_phenotype)

mat[is.na(mat)] <- 0

mat %>%
  fwrite("results/amr_out/amr_matrices/resfinder.class.decontam.csv")
