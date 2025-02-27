rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(ggpubr)
require(viridis)
require(Hmisc)

meta <- fread("data/metadata/parsed_patient_metadata.filt.csv") %>%
  filter(hap_vap_cap %in% c("HAP", "VAP"))

class_meta <- fread("results/amr_out/resfinder_out/resfinder.parsed.decontam.csv") %>%
  distinct(antimicrobial, class) %>%
  mutate(antimicrobial = capitalize(antimicrobial)) %>%
  mutate(class = capitalize(class))

genome_meta <- fread("results/assembly_out/assembly_merged_metadata.csv") %>%
  filter(gtdbtk_warnings == "N/A") %>%
  filter(run_id %in% meta$run_id) %>%
  filter(checkm2_contamination <= 5) %>%
  filter(species != "")

file_dir <- "results/amr_out/resfinder_out.assemblies/all_reports/"
file_list <- list.files(file_dir, full.names = T)

morsels <- foreach(file_name = file_list) %do% {
  id <- gsub(file_dir, "", file_name)
  id <- gsub(".txt", "", id)
  fread(file_name) %>%
    mutate(across(everything(), as.character)) %>%
    mutate(genome_name = id)
}

merged <- bind_rows(morsels) %>%
  select(genome_name, 
         antimicrobial = `# Antimicrobial`, 
         class = Class, 
         wgs_predicted_phenotype = `WGS-predicted phenotype`) %>%
  mutate(presence = wgs_predicted_phenotype == "Resistant") %>%
  inner_join(genome_meta) %>%
  mutate(antimicrobial = capitalize(antimicrobial),
         class = capitalize(class))

genome_counts <- merged %>%
  group_by(species) %>%
  summarise(n_genomes = n_distinct(genome_name)) %>%
  arrange(desc(n_genomes)) %>%
  filter(n_genomes >= 10)

amr_counts <- merged %>%
  filter(species %in% genome_counts$species) %>%
  group_by(antimicrobial, class) %>%
  summarise(n = sum(presence)) %>%
  arrange(desc(n))

# Species with most antimicrobials
merged %>% 
  distinct(checkm2_completeness, species, genome_name) %>%
  ggplot(aes(x = species, y = checkm2_completeness)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3))

summed_contig_len <- merged %>%
  distinct(genome_name, genome_size, species) %>%
  group_by(species) %>%
  summarise(sum_length = sum(genome_size))

merged %>%
  filter(wgs_predicted_phenotype == "Resistant") %>%
  group_by(species) %>%
  summarise(n_amr = n_distinct(antimicrobial)) %>%
  left_join(summed_contig_len) %>% 
  mutate(ratio = n_amr / sum_length) %>%
  left_join(genome_counts) %>%
  filter(species %in% genome_counts$species) %>%
  arrange(desc(ratio)) %>% View()


merged %>%
  filter(species %in% genome_counts$species) %>%
  mutate(species = factor(species, unique(genome_counts$species))) %>%
  mutate(antimicrobial = factor(antimicrobial, unique(amr_counts$antimicrobial))) %>%
  mutate(class = factor(class, unique(amr_counts$class))) %>%
  group_by(species, antimicrobial, class) %>%
  summarise(n = sum(presence)) %>%
  left_join(genome_counts) %>%
  mutate(prop = n / n_genomes) %>%
  filter(n > 0) %>%
  ggplot(aes(y = species, x = antimicrobial, fill = prop)) +
  geom_tile() +
  facet_grid(cols = vars(class),
             space = "free", scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3))




merged %>%
  filter(completeness)
  distinct(genome_name)
  group_by(genome_name) 
  filter(complete) %>%
  distinct(genome_name)

merged %>% 
  group_by(species, antimicrobial) %>%
  summarise(n = sum(wgs_predicted_phenotype == "Resistant")) %>%
  arrange(desc(n)) %>%
  ggplot(aes(x = species, y = antimicrobial, fill = n)) +
  geom_tile()

genome_plot_df <- merged %>%
  filter(checkm2_completeness >= 90) %>%
  group_by(genome_name) %>%
  summarise(n = sum(wgs_predicted_phenotype == "Resistant"))
  
median_resist <- median(genome_plot_df$n)

genome_plot_df %>%
  ggplot(aes(x = n)) +
  geom_histogram(fill = "#B0CBE7FF",
                 color = "black") +
  labs(x = "No. resistance types", y = "") +
  geom_vline(xintercept = median_resist, lty = "dashed")
  
  
genome_plot_df %>%
  summarise(n_multi = sum(n > 2))
merged %>%
  filter(genome_name == "INHALE_FRESH_47-barcode10-bins.2") %>% View()
merged %>%
  filt
merged %>% 
  group_by(species, `Resistance gene`, Phenotype) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>% View()


merged %>%
  filter(species == "Staphylococcus aureus")

merged %>%
  group_by(species) %>%
  summarise(n = n_distinct(run_id)) %>%
  arrange(desc(n))

merged %>%
  group_by(species, class) %>%
  summarise(n = sum(wgs_predicted_phenotype == "Resistant")) %>%
  arrange(desc(n))

merged %>% distinct(species)
merged %>% filter(species == "")
