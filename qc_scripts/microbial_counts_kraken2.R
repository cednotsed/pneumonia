setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(doParallel)

meta <- fread("data/metadata/parsed_patient_metadata.csv")
high_reads <- fread("results/qc_out/high_read_samples.csv")

rank <- "S"

dat <- fread(str_glue("results/metagenomic_out/abundance_matrix.{rank}.tsv")) %>%
  select(-any_of(c("Homo sapiens", "Homo")), -unclassified) %>%
  filter(!grepl("EXTRACTION_TEST", run)) %>%
  as_tibble() %>%
  mutate(run = gsub("INHALE_FRESH_|barcode0|barcode", "", run)) %>%
  mutate(run = gsub("a|A", "", run, ignore.case = T)) %>%
  mutate(run_id = gsub("-", "_", run)) %>%
  select(-run) %>%
  # Remove runs not in meta
  filter(run_id %in% meta$run_id) %>%
  filter(run_id %in% high_reads$run_id) %>%
  column_to_rownames("run_id") 

nrow(dat)

# Microbial counts
microbe_df <- tibble(run_id = rownames(dat),
       microbial_reads = rowSums(dat))

microbe_df %>%
  ggplot(aes(x = log10(microbial_reads))) +
  geom_histogram() +
  geom_vline(xintercept = 3,
             lty = "dashed",
             color = "red") +
  labs(x = "Log10(microbial reads)", y = "No. samples") 

microbe_df %>%
  filter(microbial_reads >= 1000) %>%
  nrow()

microbe_df %>%
  filter(microbial_reads >= 1000) %>%
  fwrite(str_glue("results/qc_out/high_microbe_samples.{rank}.csv"))
