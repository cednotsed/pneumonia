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
  distinct(run_id, bug_genus, method) %>%
  filter(bug_genus != "Negative")

bug_list <- unique(bug_df$bug_genus)
bug_list <- bug_list[bug_list %in% colnames(RA_filt)]
bug_list <- bug_list[!grepl("virus", bug_list)]

morsels <- foreach(bug_name = bug_list) %do% {
  # bug_name <- bug_list[1]
  positive_patients <- bug_df %>%
    filter(bug_genus == bug_name) %>%
    group_by(run_id) %>%
    summarise(is_cultured = sum(method == "culture") > 0)
    
  temp <- RA_filt %>%
    select(all_of(c("run_id", bug_name))) %>%
    pivot_longer(!c("run_id"), names_to = "bug", values_to = "rel_a") %>%
    inner_join(positive_patients)
  
  # temp %>% ggplot(aes(x = is_cultured, y = rel_a)) +
  #   geom_boxplot()
}


plot_df <- bind_rows(morsels)

count_df <- plot_df %>%
  group_by(is_cultured) %>%
  summarise(n = n())

plot_df %>%
  ggplot(aes(x = is_cultured, y = rel_a, fill = is_cultured)) +
  geom_boxplot() +
  geom_pwc() +
  scale_fill_manual(values = c("steelblue3", "olivedrab")) +
  theme_classic() +
  labs(x = "PCR/Culture", y = "Genus relative abundance") +
  theme(legend.position = "none")

ggsave("results/benchmarking_out/culture_positive_abundance.pdf", dpi = 600, width = 3, height = 5)

bug_df %>%
  group_by(run_id) %>%
  summarise(n_culture = sum(method == "culture"),
            n_biofire = sum(method == "biofire"),
            n_curetis = sum(method == "curetis")) %>%
filter(n_curetis > n_biofire) %>%
  nrow()

bug_df %>%
  distinct(run_id)

27/212
