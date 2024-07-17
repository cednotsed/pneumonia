rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)

meta_filt <- fread("data/metadata/parsed_patient_metadata.filt.csv") %>%
  filter(high_microbe_count) %>%
  filter(hap_vap_cap %in% c("HAP", "VAP"))

micro_meta <- fread("data/metadata/parsed_microbiology_results.csv")
df_filt <- fread("results/tax_classification_out/abundance_matrices/RA.S.zeroed.decontam.2.csv") %>%
  filter(run_id %in% meta_filt$run_id)

long_df <- df_filt %>%
  pivot_longer(!run_id, names_to = "taxa", values_to = "rel_a") %>%
  filter(rel_a != 0)

all_df <- long_df %>%
  group_by(run_id) %>%
  summarise(n_microbes = n_distinct(taxa))

all_df %>%
  ggplot(aes(x = n_microbes)) +
  geom_histogram(fill = "wheat3", color = "black") +
  theme_classic() +
  geom_vline(xintercept = median(all_df$n_microbes),
             lty = "dashed") +
  annotate(x = 10, y = 35, 
           label = str_glue("median={median(all_df$n_microbes)}"),
           geom = "text") +
  labs(x = "No. microbes identified", y = "No. patients") 

ggsave("results/benchmarking_out/num_microbes_per_patient.sequencing.pdf", width = 3, height = 3)

microbe_counts <- long_df %>%
  mutate(is_pathogen = taxa %in% c(micro_meta$bugs)) %>%
  filter(is_pathogen) %>%
  group_by(run_id) %>%
  summarise(n_microbes = n_distinct(taxa)) %>%
  right_join(meta_filt %>% select(run_id)) %>%
  mutate(n_microbes = replace_na(n_microbes, 0))

path_df <- microbe_counts %>%  
  group_by(n_microbes) %>%
  summarise(n = n_distinct(run_id))
  
path_df %>%
  mutate(n_microbes = factor(n_microbes, c(1, 2, 3, 4, 5, 0))) %>%
  mutate(prop = n / sum(path_df$n)) %>%
  ggplot(aes(x = factor(n_microbes), y = prop)) +
    geom_bar(stat = "identity",
             color = "black") +
  theme_classic() +
  labs(x = "No. pathogens identified", y = "Prop. patients") 

path_df %>%
  mutate(prop = n / sum(path_df$n)) %>%
  mutate(cat = case_when(n_microbes == 0 ~ "None", 
                         n_microbes == 1 ~ "Monomicrobial",
                         TRUE ~ "Polymicrobial")) %>%
  group_by(cat) %>%
  summarise(sum(prop),
            sum(n))
ggsave("results/benchmarking_out/polymicrobial.sequencing.pdf", width = 3, height = 3)

