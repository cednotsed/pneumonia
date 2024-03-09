rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)

meta <- fread("data/metadata/parsed_patient_metadata.csv")

high_microbe <- fread("results/qc_out/high_microbe_samples.S.csv")
high_read <- fread("results/qc_out/high_read_samples.csv")
read_counts <- fread("results/qc_out/sample_read_counts.csv")

meta_filt <- meta %>%
  filter(run_id %in% high_microbe$run_id) %>%
  filter(run_id %in% high_read$run_id) %>%
  left_join(read_counts)

dups <- meta_filt %>%
  filter(duplicated(sample_id) & 
           !grepl("water", sample_id, ignore.case = T))

dup_df <- meta_filt %>%
  filter(sample_id %in% dups$sample_id)

non_dup_df <- meta_filt %>%
  filter(!(sample_id %in% dups$sample_id))

dup_filt <- foreach(sample_name = unique(dup_df$sample_id), .combine = "bind_rows") %do% {
  dup_df %>%
    filter(sample_id == sample_name) %>%
    arrange(desc(n_reads)) %>%
    head(1)
}

final <- bind_rows(non_dup_df, dup_filt)

final %>%
  fwrite("data/metadata/parsed_patient_metadata.filt.csv")

# After filtering
final %>%
  ggplot(aes(x = log10(n_reads))) +
  geom_histogram()
