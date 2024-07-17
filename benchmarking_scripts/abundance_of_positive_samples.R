setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(doParallel)
require(ggpubr)

meta <- fread("data/metadata/parsed_patient_metadata.filt.csv") %>%
  filter(high_microbe_count)

RA_filt <- fread("results/tax_classification_out/abundance_matrices/RA.G.zeroed.decontam.2.csv") %>%
  filter(run_id %in% meta$run_id)

# Get bug list
meta <- fread("data/metadata/parsed_microbiology_results.csv")

bug_df <- meta %>%
  distinct(run_id, bug_genus) %>%
  filter(bug_genus != "Negative")

bug_list <- unique(bug_df$bug_genus)
bug_list <- bug_list[bug_list %in% colnames(RA_filt)]

morsels <- foreach(bug_name = bug_list) %do% {
  # bug_name <- bug_list[1]
  positive_patients <- bug_df %>%
    filter(bug_genus == bug_name) %>%
    distinct(run_id)
  
  temp <- RA_filt %>%
    select(all_of(c("run_id", bug_name))) %>%
    mutate(sample_type = case_when(run_id %in% positive_patients$run_id ~ "Positive",
                                   TRUE ~ "Negative")) %>%
    mutate(sample_type = factor(sample_type, c("Negative", "Positive", "Blank"))) %>%
    pivot_longer(!c("run_id", "sample_type"), names_to = "bug", values_to = "rel_a")
  
  # temp %>% ggplot(aes(x = sample_type, y = rel_a)) +
  #   geom_boxplot()
}

bind_rows(morsels) %>%
  filter(rel_a != 0) %>%
  ggplot(aes(x = sample_type, y = rel_a, fill = sample_type)) +
  geom_boxplot() +
  geom_pwc() +
  scale_fill_manual(values = c("steelblue3", "olivedrab")) +
  theme_classic() +
  labs(x = "PCR/Culture", y = "Log10(read count)") +
  theme(legend.position = "none")

bind_rows(morsels) %>%
  filter(rel_a != 0) %>%
  group_by(sample_type) %>%
  summarise(median = median(rel_a))


