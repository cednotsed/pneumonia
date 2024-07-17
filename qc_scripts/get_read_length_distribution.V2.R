rm(list = ls())
setwd("C:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)

meta <- fread("data/metadata/parsed_patient_metadata.filt.csv")

df <- fread("results/qc_out/read_lengths.txt", header = F) 

colnames(df) <- c("read_length", "path")

relevant_ids <- df %>% 
  distinct(path) %>%
  mutate(run = gsub(".fastq.gz", "", path)) %>%
  mutate(run = gsub("\\.no_human|Library|_raw|INHALE_FRESH_|barcode0|barcode", "", run)) %>%
  mutate(run = gsub("a|A", "", run, ignore.case = T)) %>%
  mutate(run_id = gsub("-", "_", run)) %>%
  filter(run_id %in% meta$run_id)

df_filt <- df %>%
  filter(path %in% relevant_ids$path) 

df_filt %>%
  summarise(median = median(read_length),
            lower = quantile(read_length, 0.25),
            upper = quantile(read_length, 0.75))

df_filt
