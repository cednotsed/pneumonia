rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)

meta_filt <- fread("data/metadata/parsed_patient_metadata.filt.csv") %>%
  filter(high_microbe_count) %>%
  filter(hap_vap_cap %in% c("HAP", "VAP"))

sp_list <- fread("data/metadata/parsed_microbiology_results.species_list.csv", sep = ",")
micro_filt <- fread("data/metadata/parsed_microbiology_results.csv") %>%
  filter(run_id %in% meta_filt$run_id)

df_filt <- fread("results/tax_classification_out/abundance_matrices/RA.S.zeroed.decontam.2.csv") %>%
  filter(run_id %in% meta_filt$run_id)

long_df <- df_filt %>%
  pivot_longer(!run_id, names_to = "taxa", values_to = "rel_a") %>%
  filter(rel_a > 0)

# Categorise micro results
seq_positive <- long_df %>%
  filter(taxa %in% sp_list$bugs) %>%
  group_by(run_id) %>%
  summarise(n = n_distinct(taxa)) %>%
  filter(n > 0)

contams <- micro_filt %>%
  filter(bugs != "Negative") %>%
  filter(!(run_id %in% long_df$run_id)) %>%
  distinct(run_id)

micro_types <- micro_filt %>%
  mutate(contams = run_id %in% contams$run_id) %>%
  mutate(rna_virus = grepl("virus", bugs) & !grepl("adenovirus", bugs)) %>%
  mutate(negative = grepl("Invalid|Negative", bugs)) %>%
  mutate(bacterial = !grepl("Invalid|Negative", bugs) & !grepl("virus", bugs)) %>%
  mutate(sequencing_positive = run_id %in% seq_positive$run_id)

type_counts <- micro_types %>%
  group_by(run_id, sequencing_positive) %>%
  summarise(type = case_when(all(contams) ~ "Only contaminants",
                             !all(contams) & any(negative) & !any(bacterial) & !any(rna_virus) ~ "PCR negative",
                             !all(contams) & any(rna_virus) & !any(bacterial) ~ "RNA virus",
                             !all(contams) & any(bacterial) ~ "Bacterial")) %>%
  mutate(sequencing_positive = run_id %in% seq_positive$run_id)

type_counts %>%
  group_by(type, sequencing_positive) %>%
  summarise(n = n_distinct(run_id))
  # filter(!sequencing_positive)

undiagnosed <- type_counts %>%
  filter(type == "PCR negative")

undiagnosed %>%
  select(run_id) %>%
  fwrite("results/benchmarking_out/undiagnosed_ids.txt")

micro_types %>%
  filter(!rna_virus, !contams) %>%
  distinct(run_id) %>%
  fwrite("results/benchmarking_out/to_assess_ids.txt")

type_counts %>%
  ungroup() %>%
  filter(!(type %in% c("RNA virus", "PCR negative"))) %>%
  summarise(n_positive = sum(sequencing_positive))

long_df

long_rel_a <- fread("results/tax_classification_out/abundance_matrices/RA.S.zeroed.decontam.2.csv") %>%
  filter(run_id %in% meta_filt$run_id) %>%
  pivot_longer(!run_id, names_to = "taxa", values_to = "rel_a") %>%
  filter(rel_a > 0)

long_read <- fread("results/tax_classification_out/abundance_matrices/read_counts.S.zeroed.decontam.2.csv") %>%
  filter(run_id %in% meta_filt$run_id) %>%
  pivot_longer(!run_id, names_to = "taxa", values_to = "read_count") %>%
  filter(read_count > 0)

long_df <- df_filt %>%
  pivot_longer(!run_id, names_to = "taxa", values_to = "rel_a") %>%
  filter(rel_a > 0)
test <- type_counts %>%
  filter(type == "PCR negative") %>%
  filter(sequencing_positive)

long_read %>%
  filter(run_id %in% test$run_id) %>%
  filter(taxa == "Escherichia coli") %>%
  left_join(long_rel_a %>% select(run_id, taxa, rel_a))


# # Get dominant microbes
# undiagnosed 
# long_df %>%
#   filter(run_id)
# morsels <- foreach(id = unique(long_df$run_id)) %do% {
#   temp <- long_df %>%
#     filter(run_id == id) %>%
#     arrange(desc(rel_a))
#   
#   temp %>%
#     mutate(dominant = c(T, rep(F, nrow(temp) - 1)))
# }
# 
# parsed <- bind_rows(morsels)
# # filter(!grepl("Candida|Aspergillus|Nakaseomyces", taxa)) 
# 
# count_df <- parsed %>%
#   group_by(taxa) %>%
#   summarise(n = n_distinct(run_id)) %>%
#   arrange(desc(n)) %>%
#   mutate(prop = n / nrow(meta_filt)) %>%
#   head(30)
# 
# plot_df <- parsed %>%
#   filter(taxa %in% count_df$taxa) %>%
#   group_by(taxa) %>%
#   summarise(n = n_distinct(run_id)) %>%
#   arrange(desc(n)) %>%
#   mutate(prop = n / nrow(meta_filt))
# 
# plot_df %>%  
#   mutate(taxa = factor(taxa, count_df$taxa)) %>%
#   ggplot(aes(x = taxa, y = prop)) +
#   geom_bar(stat = "identity", color = "black") +
#   theme_classic() +
#   scale_fill_manual(values = c("darkseagreen4", "khaki3")) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3,
#                                    face = "italic")) +
#   labs(x = "Microbial species", y = "Prop. patients")
# 
# ggsave("results/benchmarking_out/bug_proportions.no_pathogens.pdf", 
#        dpi = 600, 
#        width = 7, 
#        height = 3)
# 
# virus <- micro_filt %>%
#   filter(!(bugs %in% c("Negative", "Invalid test"))) %>%
#   filter(grepl("virus", bugs)) %>% 
#   filter(!grepl("adenovirus", bugs)) %>%
#   distinct(run_id)
# 
# pcr_positive <- micro_filt %>%
#   filter(!(bugs %in% c("Negative", "Invalid test"))) %>%
#   filter(!(run_id %in% virus$run_id)) %>%
#   distinct(run_id)
# 
# negatives <- meta_filt %>%
#   filter(!(run_id %in% c(virus$run_id, pcr_positive$run_id))) %>%
#   distinct(run_id)
# 
# 
# neg_plot_df <- long_df %>% 
#   filter(run_id %in% negatives$run_id) %>% 
#   separate(taxa, c("genus"), "\\ ") %>%
#   group_by(run_id, genus) %>%
#   summarise(rel_a = sum(rel_a)) %>%
#   filter(rel_a > 0.1)
# 
# pal <- distinctColorPalette(n_distinct(neg_plot_df$genus))
# 
# neg_plot_df %>%  
#   ggplot(aes(x = run_id, y = rel_a, fill = genus)) +
#   geom_bar(stat = "identity", color = "black") +
#   scale_fill_manual(values = pal)
