rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)

meta_filt <- fread("data/metadata/parsed_patient_metadata.filt.csv") %>%
  filter(high_microbe_count) %>%
  filter(hap_vap_cap %in% c("HAP", "VAP"))

df_filt <- fread("results/tax_classification_out/abundance_matrices/RA.S.zeroed.decontam.2.csv") %>%
  filter(run_id %in% meta_filt$run_id)

long_df <- df_filt %>%
  pivot_longer(!run_id, names_to = "taxa", values_to = "rel_a") %>%
  filter(rel_a != 0)

plot_df <- long_df %>%
  group_by(taxa) %>%
  summarise(n = n_distinct(run_id)) %>%
  arrange(desc(n)) %>%
  head(30)

plot_df %>%  
  mutate(taxa = factor(taxa, plot_df$taxa)) %>%
  ggplot(aes(x = taxa, y = n)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Get dominant pathogen
morsels <- foreach(id = unique(long_df$run_id)) %do% {
  long_df %>%
    filter(run_id == id) %>%
    arrange(desc(rel_a)) %>%
    head(1)
}

merged <- bind_rows(morsels) %>%
  mutate(above_half = rel_a > 0.5)

count_df <- merged %>%
  group_by(taxa) %>%
  summarise(n = n_distinct(run_id)) %>%
  filter(!(taxa %in% c("Candida albicans", "Aspergillus fumigatus", "Nakaseomyces glabratus"))) %>%
  arrange(desc(n)) %>%
  head(30)

plot_df <- merged %>%
  filter(taxa %in% count_df$taxa) %>%
  group_by(taxa, above_half) %>%
  summarise(n = n_distinct(run_id))

plot_df %>%
  filter(taxa == "Escherichia coli")

plot_df %>%  
  mutate(prop = n / sum(plot_df$n)) %>%
  mutate(taxa = factor(taxa, unique(count_df$taxa))) %>%
  ggplot(aes(x = taxa, y = n, fill = above_half)) +
  geom_bar(color = "black", stat = "identity") +
  theme_classic() +
  scale_fill_manual(values = c("olivedrab", "goldenrod")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Bacterial species", y = "Prop. patients")

merged %>%
  filter(taxa %in% count_df$taxa) %>%
  mutate(taxa = factor(taxa, unique(count_df$taxa))) %>%
  mutate(is_commensal = !(taxa %in% c("Escherichia coli", "Staphylococcus aureus",
                                   "Elizabethkingia miricola", "Pseudomonas aeruginosa",
                                   "Klebsiella pneumoniae", "Haemophilus influenzae",
                                   "Paraburkholderia fungorum", "Citrobacter koseri",
                                   "Enterobacter hormaechei", "Enterococcus faecium",
                                   "Moraxella catarrhalis", "Proteus mirabilis",
                                   "Acinetobacter baumannii", "Morganella morganii",
                                   "Stenotrophomonas maltophilia", "Citrobacter freundii",
                                   "Klebsiella aerogenes", "Achromobacter xylosoxidans"))) %>%
  ggplot(aes(x = taxa, y = rel_a, fill = is_commensal)) +
  geom_boxplot() +
  facet_grid(cols = vars(is_commensal),
             space = "free", scales = "free") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  
  
bind_rows(morsels) %>%
  filter(taxa %in% plot_df$taxa) %>%
  filter(!(taxa %in% c("Candida", "Aspergillus", "Nakaseomyces"))) %>%
  mutate(taxa = factor(taxa, unique(plot_df$taxa))) %>%
  ggplot(aes(x = taxa, y = rel_a)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

long_df %>%
  filter(taxa == "Elizabethkingia") %>% View()
