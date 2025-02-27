rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(ggpubr)
require(viridis)
require(Hmisc)

meta <- fread("data/metadata/parsed_patient_metadata.filt.csv") %>%
  filter(hap_vap_cap %in% c("HAP", "VAP"))

class_meta <- fread("results/amr_out/resfinder_out/resfinder.parsed.decontam.csv") %>%
  distinct(antimicrobial, class) %>%
  mutate(antimicrobial = capitalize(antimicrobial)) %>%
  mutate(class = capitalize(class))

genome_meta <- fread("results/assembly_out/assembly_merged_metadata.csv") %>%
  filter(gtdbtk_warnings == "N/A") %>%
  filter(species != "") %>%
  filter(run_id %in% meta$run_id) %>%
  filter(checkm2_contamination <= 5) %>%
  mutate(complete = checkm2_completeness >= 90)

genome_meta %>%
  mutate(complete = checkm2_completeness >= 90) %>%
  group_by(complete) %>%
  summarise(n = n())

file_dir <- "results/amr_out/resfinder_out.assemblies/all_tables/"
file_list <- list.files(file_dir, full.names = T)

morsels <- foreach(file_name = file_list) %do% {
  id <- gsub(file_dir, "", file_name)
  id <- gsub(".txt", "", id)
  fread(file_name) %>%
    mutate(across(everything(), as.character)) %>%
    mutate(genome_name = id)
}

merged <- bind_rows(morsels) %>%
  inner_join(genome_meta) %>%
  filter(run_id %in% meta$run_id) %>%
  filter(complete)

genome_n

genome_meta %>%
  filter(complete) %>%
  filter(!(genome_name %in% merged$genome_name))
  nrow()
  group_by(genome_name) %>%
  sum
