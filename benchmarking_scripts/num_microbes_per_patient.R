rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)

meta_filt <- fread("data/metadata/parsed_patient_metadata.filt.csv") %>%
  filter(hap_vap_cap %in% c("HAP", "VAP"))

micro_meta <- fread("data/metadata/parsed_microbiology_results.csv", sep = ",") %>%
  filter(run_id %in% meta_filt$run_id)

total <- n_distinct(micro_meta$run_id)

parsed <- micro_meta %>%
  filter(!(bug_genus %in% c("Invalid", "Negative"))) %>%
  distinct(run_id, bug_genus) 
  
plot_df <- parsed %>%
  group_by(run_id) %>%
  summarise(n_genus = n_distinct(bug_genus)) %>%
  right_join(micro_meta %>% distinct(run_id)) %>% 
  mutate(n_genus = replace_na(n_genus, 0)) %>%
  group_by(n_genus) %>%
  summarise(n = n_distinct(run_id)) %>%
  mutate(prop = n / total)

plot_df %>%
  mutate(n_genus = factor(n_genus, c(1, 2, 3, 4, 5, 0))) %>%
  ggplot(aes(x = n_genus, y = prop)) +
  geom_col(color = "black") + 
  theme_classic() +
  labs(x = "No. pneumonia-associated genera", y = "Prop. patients")

ggsave("results/benchmarking_out/polymicrobial.pdf",
       dpi = 600, 
       width = 7, 
       height = 3)
