rm(list = ls())
setwd("C:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)

micro_meta <- fread("data/metadata/parsed_microbiology_results.csv") %>%
  filter(!(bugs %in% c("Invalid test"))) %>%
  mutate(bugs = ifelse(grepl("Enterobacter", bugs), "Enterobacter spp.", bugs)) %>%
  mutate(bugs = ifelse(grepl("Acinetobacter", bugs), "Acinetobacter spp.", bugs)) %>%
  mutate(bugs = ifelse(grepl("Proteus", bugs), "Proteus spp.", bugs)) %>%
  filter(bugs != "Coliform")

meta_filt <- fread("data/metadata/parsed_patient_metadata.filt.csv") %>%
  filter(hap_vap_cap %in% c("HAP", "VAP")) %>%
  filter(high_microbe_count)

# Parse sequencing results
control_meta <- fread("data/metadata/parsed_patient_metadata.filt.csv") %>%
  filter(hap_vap_cap == "Water control")

RA <- fread("results/tax_classification_out/abundance_matrices/RA.S.zeroed.csv") 

RA_filt <- RA %>%
  filter(run_id %in% meta_filt$run_id)

control_RA <- RA %>%
  filter(run_id %in% control_meta$run_id)

read <- fread("results/tax_classification_out/abundance_matrices/read_counts.S.zeroed.csv")

read_filt <- read %>%
  filter(run_id %in% meta_filt$run_id)

control_read <- read %>%
  filter(run_id %in% control_meta$run_id)

long_df <- RA_filt %>%
  pivot_longer(!run_id, names_to = "species", values_to = "rel_a") %>%
  left_join(meta_filt) %>%
  filter(high_microbe_count)

long_read_df <- read_filt %>%
  pivot_longer(!run_id, names_to = "species", values_to = "read_count") %>%
  left_join(meta_filt) %>%
  filter(high_microbe_count)

control_read_df <- control_read %>%
  pivot_longer(!run_id, names_to = "species", values_to = "read_count") %>%
  left_join(control_meta) %>%
  filter(high_microbe_count) %>%
  select(species, run, control_reads = read_count)

control_RA_df <- control_RA %>%
  pivot_longer(!run_id, names_to = "species", values_to = "rel_a") %>% 
  left_join(control_meta) %>%
  select(species, run, control_RA = rel_a) %>%
  group_by(species, run) %>%
  summarise(control_RA = max(control_RA))

sequencing_df <- long_df %>%
  arrange(desc(species)) %>%
  left_join(long_read_df %>% select(run_id, species, read_count)) %>%
  left_join(control_RA_df) %>%
  # Remove pathogens that are found in water controls
  filter(control_RA < 5 * rel_a)

sequencing_df

# Culture and PCR 
total_counts <- micro_meta %>%
  group_by(method) %>%
  summarise(n_total = n_distinct(run_id))

# plot_df <- micro_meta %>%
#   group_by(method, bugs) %>%
#   summarise(n = n_distinct(run_id)) %>%
#   left_join(total_counts) %>%
#   mutate(prop = n / n_total * 100) %>%
#   mutate(pathogen_type = ifelse(grepl("virus", bugs), "Viral", "Bacterial")) %>%
#   arrange(desc(prop))

plot_df <- micro_meta %>%
  mutate(pathogen_type = ifelse(grepl("virus", bugs), "Viral", "Bacterial")) %>%
  mutate(bugs = gsub("Escherichia", "E.", bugs)) %>%
  mutate(bugs = gsub("Staphylococcus", "S.", bugs)) %>%
  mutate(bugs = gsub("Pseudomonas", "P.", bugs)) %>%
  mutate(bugs = gsub("Human rhinovirus/Enterovirus", "Rhinovirus/enterovirus", bugs)) %>%
  mutate(bugs = gsub("Klebsiella", "K.", bugs)) %>%
  mutate(bugs = gsub("Haemophilus", "H.", bugs)) %>%
  mutate(bugs = gsub("Alphainfluenzavirus influenzae", "Influenza A", bugs)) %>%
  mutate(bugs = gsub("Betainfluenzavirus influenzae", "Influenza B", bugs)) %>%
  mutate(bugs = gsub("Stenotrophomonas", "Steno.", bugs)) %>%
  mutate(bugs = gsub("Streptococcus", "Strep.", bugs)) %>%
  mutate(bugs = gsub("Human metapneumo", "Metapneumo", bugs)) %>%
  mutate(bugs = gsub("Respiratory syncytial virus", "RSV", bugs)) %>%
  group_by(bugs, pathogen_type) %>%
  filter(bugs != "Negative") %>%
  summarise(n = n_distinct(run_id)) %>%
  mutate(prop = n / n_distinct(micro_meta$run_id) * 100) %>%
  arrange(desc(n))

plot_df %>%  
  mutate(bugs = factor(bugs, unique(plot_df$bugs))) %>%
  ggplot(aes(x = bugs, y = prop, fill = pathogen_type)) +
  geom_bar(stat = "identity", color = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3,
                                   face = "italic")) +
  labs(x = "Pathogen", y = "% of positive tests", fill = "Pathogen type") 
  
ggsave("results/benchmarking_out/bug_proportions.pdf", 
       dpi = 600, 
       width = 7, 
       height = 3)
  
  
