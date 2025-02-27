rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(doParallel)

meta <- fread("data/metadata/metadata_221024.csv") %>%
  rename_all(~tolower(gsub(" |/", "_", .x))) %>%
  mutate(run = as.character(run)) %>%
  mutate(biofire_organism = ifelse(biofire_validity == "Valid", biofire_organism, "Invalid test")) %>%
  mutate(curetis_organism = ifelse(curetis_validity == "Valid", curetis_organism, "Invalid test")) %>%
  mutate(culture_organism = gsub("\n", ";", culture_organism)) %>%
  mutate(biofire_organism = gsub("\n", ";", biofire_organism)) %>%
  mutate(curetis_organism = gsub("\n", ";", curetis_organism)) %>%
  # Parse viruses
  mutate(culture_organism = gsub("Adenovirus", "Human Adenovirus", culture_organism)) %>%
  mutate(culture_organism = gsub("Influenza A", "Alphainfluenzavirus influenzae", culture_organism)) %>%
  mutate(culture_organism = gsub("Influenza B", "Betainfluenzavirus influenzae", culture_organism)) %>%
  mutate(biofire_organism = gsub("Adenovirus", "Human Adenovirus", biofire_organism)) %>%
  mutate(biofire_organism = gsub("Influenza A", "Alphainfluenzavirus influenzae", biofire_organism)) %>%
  mutate(biofire_organism = gsub("Influenza B", "Betainfluenzavirus influenzae", biofire_organism)) %>%
  mutate(curetis_organism = gsub("Adenovirus", "Human Adenovirus", curetis_organism)) %>%
  mutate(curetis_organism = gsub("Influenza A", "Alphainfluenzavirus influenzae", curetis_organism)) %>%
  mutate(curetis_organism = gsub("Influenza B", "Betainfluenzavirus influenzae", curetis_organism)) %>%
  mutate(culture_organism = gsub("Enterobacter aerogenes", "Klebsiella aerogenes", culture_organism)) %>%
  mutate(curetis_organism = gsub("Enterobacter aerogenes", "Klebsiella aerogenes", curetis_organism)) %>%
  mutate(biofire_organism = gsub("Enterobacter aerogenes", "Klebsiella aerogenes", biofire_organism)) %>%
  mutate(culture_organism = gsub("Raoultella ornitholytica", "Raoultella ornithinolytica", culture_organism)) %>%
  mutate(curetis_organism = gsub("Raoultella ornitholytica", "Raoultella ornithinolytica", curetis_organism)) %>%
  mutate(biofire_organism = gsub("Raoultella ornitholytica", "Raoultella ornithinolytica", biofire_organism)) %>%
  mutate(culture_organism = gsub("Streptococcus pnuemoniae", "Streptococcus pneumoniae", culture_organism)) %>%
  mutate(curetis_organism = gsub("Streptococcus pnuemoniae", "Streptococcus pneumoniae", curetis_organism)) %>%
  mutate(biofire_organism = gsub("Streptococcus pnuemoniae", "Streptococcus pneumoniae", biofire_organism)) %>%
  mutate(culture_organism = gsub("Haemophilus influenzae", "Haemophilus influenza", culture_organism)) %>%
  mutate(curetis_organism = gsub("Haemophilus influenzae", "Haemophilus influenza", curetis_organism)) %>%
  mutate(biofire_organism = gsub("Haemophilus influenzae", "Haemophilus influenza", biofire_organism)) %>%
  mutate(culture_organism = gsub("Haemophilus influenza", "Haemophilus influenzae", culture_organism)) %>%
  mutate(curetis_organism = gsub("Haemophilus influenza", "Haemophilus influenzae", curetis_organism)) %>%
  mutate(biofire_organism = gsub("Haemophilus influenza", "Haemophilus influenzae", biofire_organism)) %>%
  mutate(run_id = str_glue("{run}_{barcode}")) %>%
  mutate(hap_vap_cap = ifelse(grepl("Water", sample_id, ignore.case = T),
                              "Water control",
                              hap_vap_cap)) %>%
  filter(!(sample_id %in% c("Special", "Failed Run - 26"))) %>%
  mutate(run_id = gsub("a|A", "", run_id))

# Controls
control <- fread("data/metadata/control_metadata.csv") %>%
  mutate(barcode = as.character(barcode)) %>%
  mutate(run_id = str_glue("{run}_{barcode}")) %>%
  mutate(hap_vap_cap = "Healthy", hospital = "Local GP (controls)") %>%
  mutate(sample_type = "SPU") %>%
  mutate(hap_vap_cap = ifelse(run == "A", "Water control", hap_vap_cap),
         hospital = ifelse(run == "A", NA, hospital))

meta %>%
  bind_rows(control) %>% 
  mutate(hap_vap2 = case_when(hap_vap_cap == "HAP" & vent_length_hours <= 0 ~ "NV-HAP",
                              hap_vap_cap == "HAP" & vent_length_hours > 0 ~ "V-HAP",
                              hap_vap_cap == "VAP" ~ "VAP",
                              hap_vap_cap == "healthy" ~ "Healthy",
                              hap_vap_cap == "Water control" ~ "Water control")) %>%
  mutate(ventilation = vent_length_hours > 0) %>%
  fwrite("data/metadata/parsed_patient_metadata.csv")


