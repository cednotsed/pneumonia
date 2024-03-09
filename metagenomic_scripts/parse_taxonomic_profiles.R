setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(doParallel)

## FUNCTIONS ##
otu_to_RA <- function(df) {
  mat <- as.matrix(df)
  RA_df <- as.data.frame(mat / rowSums(mat))
  # RA_df <- add_column(RA_df, df$run, .before = 1)
  colnames(RA_df) <- colnames(df)
  
  return(RA_df)
}

otu_to_PA <- function(df, read_threshold) {
  prev_read <- df
  
  prev_read[prev_read <= read_threshold] <- 0
  prev_read[prev_read > read_threshold] <- 1
  colnames(prev_read) <- colnames(df)
  return(prev_read)
}

RA_to_PA <- function(RA_df, PA_threshold) {
  prev_RA <- RA_df
  prev_RA[prev_RA <= PA_threshold] <- 0
  prev_RA[prev_RA > PA_threshold] <- 1
  colnames(prev_RA) <- colnames(RA_df)
  return(prev_RA)
}

rank <- "F"

meta <- fread("data/metadata/parsed_patient_metadata.csv")
high_microbes <- fread(str_glue("results/qc_out/high_microbe_samples.{rank}.csv"))

dat <- fread(str_glue("results/metagenomic_out/abundance_matrix.{rank}.tsv")) %>%
  select(-any_of(c("Homo sapiens", "Homo", "Hominidae")), -unclassified) %>%
  filter(!grepl("EXTRACTION_TEST", run)) %>%
  as_tibble() %>%
  mutate(run = gsub("INHALE_FRESH_|barcode0|barcode", "", run)) %>%
  mutate(run = gsub("a|A", "", run, ignore.case = T)) %>%
  mutate(run_id = gsub("-", "_", run)) %>%
  select(-run) %>%
  # Remove runs not in meta
  filter(run_id %in% meta$run_id) %>%
  filter(run_id %in% high_microbes$run_id) %>%
  column_to_rownames("run_id") 

# Filter read counts <= 10, RA <= 0.005
RA_df <- otu_to_RA(dat)
prev_RA <- RA_to_PA(RA_df, 0.005)
prev_read <- otu_to_PA(dat, 10)
prev_df <- as.data.frame(prev_read & prev_RA)

# Set read counts to zero
dat_filt <- dat
dat_filt[!as.matrix(prev_df)] <- 0

# Remove zero taxa and IDs with zero reads
dat_filt <- dat_filt[rowSums(dat_filt) != 0, colSums(dat_filt) != 0]

# Get relative abundance
RA_filt <- otu_to_RA(dat_filt)

# Save files
dat_filt %>% 
  rownames_to_column("run_id") %>%
  fwrite(str_glue("results/metagenomic_out/read_counts.{rank}.zeroed.csv"),
         row.names = F)

RA_filt %>% 
  rownames_to_column("run_id") %>%
  fwrite(str_glue("results/metagenomic_out/RA.{rank}.zeroed.csv"),
         row.names = F)

nrow(RA_filt)
# dat_filt %>%
#   filter(run_id == "1_1")
# otu_to_RA(dat_filt) %>%
#   rownames_to_column("run_id") %>%
#   filter(run_id == "1_1") %>%
#   pivot_longer(!run_id, names_to = "species", values_to = "rel_a") %>%
#   arrange(desc(rel_a))
# select(`Pseudomonas aeruginosa`)
# 
# dat_filt["1_1", ] %>%
#   pivot_longer(everything(), names_to = "species", values_to = "read_count") %>%
#   arrange(desc(read_count))
