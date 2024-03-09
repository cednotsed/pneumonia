rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)

meta <- fread("data/metadata/parsed_patient_metadata.csv")

meta %>%
  mutate(file_name = gsub("_", "-", run_id)) %>%
  mutate(file_name = gsub("-", "-barcode", file_name)) %>%
  mutate(file_name = gsub("barcode12", "barcode12a", file_name)) %>%
  mutate(file_name = ifelse(grepl("barcode10|barcode11|barcode12", file_name), 
                            file_name,
                            gsub("barcode", "barcode0", file_name))) %>%
  mutate(file_name = str_glue("INHALE_FRESH_{file_name}.fastq.gz")) %>%
  select(file_name) %>%
  fwrite("data/metadata/patient_file_paths.txt",
         col.names = F,
         row.names = F,
         eol = "\n")

meta %>%
  filter(run_id == "15_12")
