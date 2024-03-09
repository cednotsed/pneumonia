setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(doParallel)

meta <- fread("data/metadata/patient_metadata.csv") %>%
  rename_all(~tolower(gsub(" |/", "_", .x))) %>%
  mutate(run = as.character(run)) %>%
  mutate(biofire_organism = ifelse(biofire_valid == "Valid", biofire_organism, "Invalid test")) %>%
  mutate(curetis_organism = ifelse(curetis_valid == "Valid", curetis_organism, "Invalid test")) %>%
  mutate(micro_organism = gsub("\n", ";", micro_organism)) %>%
  mutate(biofire_organism = gsub("\n", ";", biofire_organism)) %>%
  mutate(curetis_organism = gsub("\n", ";", curetis_organism)) %>%
  # Parse viruses
  mutate(micro_organism = gsub("Adenovirus", "Human Adenovirus", micro_organism)) %>%
  mutate(micro_organism = gsub("Influenza A", "Alphainfluenzavirus influenzae", micro_organism)) %>%
  mutate(micro_organism = gsub("Influenza B", "Betainfluenzavirus influenzae", micro_organism)) %>%
  mutate(biofire_organism = gsub("Adenovirus", "Human Adenovirus", biofire_organism)) %>%
  mutate(biofire_organism = gsub("Influenza A", "Alphainfluenzavirus influenzae", biofire_organism)) %>%
  mutate(biofire_organism = gsub("Influenza B", "Betainfluenzavirus influenzae", biofire_organism)) %>%
  mutate(curetis_organism = gsub("Adenovirus", "Human Adenovirus", curetis_organism)) %>%
  mutate(curetis_organism = gsub("Influenza A", "Alphainfluenzavirus influenzae", curetis_organism)) %>%
  mutate(curetis_organism = gsub("Influenza B", "Betainfluenzavirus influenzae", curetis_organism)) %>%
  mutate(micro_organism = gsub("Enterobacter aerogenes", "Klebsiella aerogenes", micro_organism)) %>%
  mutate(curetis_organism = gsub("Enterobacter aerogenes", "Klebsiella aerogenes", curetis_organism)) %>%
  mutate(biofire_organism = gsub("Enterobacter aerogenes", "Klebsiella aerogenes", biofire_organism)) %>%
  mutate(micro_organism = gsub("Raoultella ornitholytica", "Raoultella ornithinolytica", micro_organism)) %>%
  mutate(curetis_organism = gsub("Raoultella ornitholytica", "Raoultella ornithinolytica", curetis_organism)) %>%
  mutate(biofire_organism = gsub("Raoultella ornitholytica", "Raoultella ornithinolytica", biofire_organism)) %>%
  mutate(micro_organism = gsub("Streptococcus pnuemoniae", "Streptococcus pneumoniae", micro_organism)) %>%
  mutate(curetis_organism = gsub("Streptococcus pnuemoniae", "Streptococcus pneumoniae", curetis_organism)) %>%
  mutate(biofire_organism = gsub("Streptococcus pnuemoniae", "Streptococcus pneumoniae", biofire_organism)) %>%
  mutate(micro_organism = gsub("Haemophilus influenzae", "Haemophilus influenza", micro_organism)) %>%
  mutate(curetis_organism = gsub("Haemophilus influenzae", "Haemophilus influenza", curetis_organism)) %>%
  mutate(biofire_organism = gsub("Haemophilus influenzae", "Haemophilus influenza", biofire_organism)) %>%
  mutate(micro_organism = gsub("Haemophilus influenza", "Haemophilus influenzae", micro_organism)) %>%
  mutate(curetis_organism = gsub("Haemophilus influenza", "Haemophilus influenzae", curetis_organism)) %>%
  mutate(biofire_organism = gsub("Haemophilus influenza", "Haemophilus influenzae", biofire_organism)) %>%
  mutate(run_id = str_glue("{run}_{barcode}")) %>%
  mutate(hap_vap_cap = ifelse(grepl("Water", sample_id, ignore.case = T),
                              "Water control",
                              hap_vap_cap)) %>%
  filter(!(sample_id %in% c("Special", "Failed Run - 26"))) %>%
  mutate(run_id = gsub("a|A", "", run_id))


fwrite(meta, "data/metadata/parsed_patient_metadata.csv")

