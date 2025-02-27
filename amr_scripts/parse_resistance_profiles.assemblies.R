rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(vegan)
require(hillR)
require(ggpubr)
require(randomcoloR)

ass_meta <- fread("results/assembly_out/assembly_merged_metadata.csv")

meta <- fread("data/metadata/parsed_patient_metadata.filt.csv") %>%
  filter(high_microbe_count)

file_dir <- "results/amr_out/resfinder_out.assemblies/all_reports/"
file_list <- list.files(file_dir, full.names = T)

res_df <- foreach(file = file_list, .combine = "bind_rows") %do% {
  id <- gsub(file_dir, "", file)
  
  fread(file) %>%
    mutate(run_name = id)
}

parsed <- res_df %>%
  rename_all(~tolower(gsub(" |-", "_", .x))) %>%
  dplyr::rename(antimicrobial = `#_antimicrobial`) %>%
  mutate(genome_name = gsub(".txt", "", run_name)) %>%
  mutate(run_id = gsub("INHALE_FRESH_|.no_human|Library|_raw", "", run_name)) %>%
  mutate(run_id = gsub(".pheno_table.txt", "", run_id)) %>%
  mutate(run_id = gsub("barcode", "", run_id)) %>%
  mutate(run_id = gsub("-0", "-", run_id)) %>%
  mutate(run_id = gsub("-", "_", run_id)) %>%
  mutate(run_id = gsub("12a", "12", run_id)) %>%
  separate(run_id, c("run", "index"), "_") %>%
  mutate(run_id = str_glue("{run}_{index}")) %>%
  inner_join(meta %>% select(run_id, hap_vap_cap, run)) %>%
  select(-run_name) %>%
  arrange(desc(antimicrobial)) %>%
  left_join(ass_meta %>% select(genome_name, family, genus, 
                                species, checkm2_completeness, checkm2_contamination, 
                                genome_size, gtdbtk_warnings))

parsed %>%
  fwrite("results/amr_out/resfinder_out/resfinder.assemblies.parsed.csv")

