setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(doParallel)

meta <- fread("data/metadata/parsed_patient_metadata.filt.csv")


rank <- "G"

dat <- fread(str_glue("results/tax_classification_out/abundance_matrices/abundance_matrix.{rank}.tsv")) %>%
  select(-any_of(c("Homo sapiens", "Homo")), -unclassified) %>%
  as_tibble() %>%
  column_to_rownames("run_id") 

rownames(dat)

# Microbial counts
microbe_df <- tibble(run_id = rownames(dat),
       microbial_reads = rowSums(dat))

# microbe_df %>% arrange(microbial_reads) %>% View()
microbe_df %>%
  ggplot(aes(x = log10(microbial_reads))) +
  geom_histogram() +
  geom_vline(xintercept = 2,
             lty = "dashed",
             color = "red") +
  labs(x = "Log10(microbial reads)", y = "No. samples") 

microbe_df %>%
  # filter(microbial_reads >= 100) %>%
  nrow()

microbe_df %>%
  # filter(microbial_reads >= 100) %>%
  fwrite(str_glue("results/qc_out/high_microbe_samples.{rank}.csv"))
