rm(list = ls())
setwd("C:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)

meta <- fread("data/metadata/parsed_patient_metadata.csv")

file_dir <- "data/sequencing_summaries/"
file_list <- list.files(file_dir, full.names = T)

morsels <- foreach(file = file_list) %do% {
  # file = file_list[1]
  id <- gsub(file_dir, "", file)
  id <- gsub("_sequencing_summary.parsed.csv.gz", "", id)
  
  temp <- fread(file) %>%
    filter(passes_filtering) %>%
    filter(barcode_arrangement != "unclassified") %>%
    select(read_id, sequence_length_template, barcode_arrangement)
  
  temp %>%
    mutate(run = id)
}

# PASS FAIL UNCLASSIFIED
read_df <- bind_rows(morsels) %>%
  mutate(run = gsub("INHALE_FRESH_", "", run)) %>%
  mutate(barcode_arrangement = gsub("barcode0|barcode", "", barcode_arrangement)) %>%
  mutate(run_id = str_glue("{run}_{barcode_arrangement}"))

read_df %>%
  mutate(run = gsub("INHALE_FRESH_|barcode0|barcode", "", run)) %>%
  mutate(run = gsub("a|A", "", run, ignore.case = T)) %>%
  mutate(run_id = gsub("-", "_", run)) %>%
  select(-run)

meta <- fread("data/metadata/parsed_patient_metadata.csv")
meta
read_df %>%
  filter(run_id %in% meta$run_id) %>%
  summarise(mean = mean(sequence_length_template),
            sd = sd(sequence_length_template))
  
read_df 
  mutate(run = gsub("INHALE_FRESH_", "", run)) %>%
  arrange(desc(total_reads))
