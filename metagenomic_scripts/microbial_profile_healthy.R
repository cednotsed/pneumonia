rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(vegan)
require(hillR)
require(ggpubr)
require(randomcoloR)

meta <- fread("data/metadata/parsed_patient_metadata.filt.csv") %>%
  filter(high_microbe_count)

controls <- meta %>%
  filter(hap_vap_cap == "Water control")

healthy <- meta %>%
  filter(hap_vap_cap == "healthy")

RA <- fread("results/tax_classification_out/abundance_matrices/RA.G.zeroed.csv")

irep_df <- fread("results/irep_out/irep_results.parsed.tsv")

irep_binary <- irep_df %>%
  separate(taxa, c("genus"), "\\ ") %>%
  group_by(run_id, genus) %>%
  summarise(replicating = any(!is.na(bPTR)))

control_RA <- RA %>%
  filter(run_id %in% controls$run_id)

healthy_RA <- RA %>%
  filter(run_id %in% healthy$run_id)

long_healthy <- healthy_RA %>%
  pivot_longer(!run_id, names_to = "taxa", values_to = "rel_a") %>%
  filter(rel_a != 0)

long_control <- control_RA %>%
  pivot_longer(!run_id, names_to = "taxa", values_to = "rel_a") %>%
  filter(rel_a != 0) %>%
  group_by(rel_a != 0)

healthy_mat <- long_healthy %>%
  # filter(!(taxa %in% long_control$taxa)) %>% View()
  pivot_wider(id_cols = run_id, names_from = taxa, values_from = rel_a) %>%
  mutate(across(everything(), ~replace_na(., 0)))

# HCLUST
healthy_mat <- healthy_mat %>% column_to_rownames("run_id")
df_dist <- vegdist(healthy_mat, method = "bray")
clustering <- hclust(df_dist)

pal <- distinctColorPalette(n_distinct(long_healthy$taxa))
pal <-  c("#DF58D0", "#DBBDD9", "#C1E482", "#DE7D52", "#DD557D", "#E2B04B", "#79BA81", "#B2EC63",
          "#65EED1", "#8D7ADD", "#DCDE49", "#9FE9D9", "#955688", "#759C9F", "#E1D6BA", "#9748E6",
          "#86EC34", "#DF9C9E", "#6FE792", "#C9F0B3", "#D8E2E6", "#82A1DB", "#878755", "#56DE58",
          "#E5D587", "#7ECEE5", "#DF9CE1")

long_healthy %>%
  left_join(meta) %>%
  mutate(run_id = factor(run_id, rownames(healthy_mat)[clustering$order])) %>%
  ggplot(aes(x = run_id, y = rel_a, fill = taxa)) +
  geom_bar(stat = "identity",
           position = "stack",
           color = "black") +
  scale_fill_manual(values = pal, na.translate = F) +
  theme_bw() +
  # facet_wrap(~hap_vap_cap, 
  #            nrow = 2, ncol = 2, scales = "free") +
  theme(axis.text.x = element_blank(),
        # axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position = "top") +
  labs(x = "Patient", y = "Relative abundance",
       fill = "Genus")

geom_text(aes(x = 0, y = Inf, 
              fill = NA, 
              label = str_glue("n={n}")),
          data = count_df,
          vjust = 1.1,
          hjust = -0.1)

ggsave("results/metagenomic_out/relative_abundance.positive_bugs.pdf", width = 12, height = 5)


foreach()
long_filt <- long_df %>%
  filter(hap_vap_cap == "VAP") %>%
  filter(rel_a > 0.5)

long_filt %>%
  group_by(genus) %>%
  summarise(n_samples = n_distinct(run_id)) %>%
  arrange(desc(n_samples))

set.seed(66)

morsels <- foreach(i = seq(1000)) %do% {
  long_filt %>%
    sample_n(60, replace = F) %>%
    group_by(genus) %>%
    summarise(n_samples = n_distinct(run_id)) %>%
    mutate(index = i)
}

n_distinct(long_filt$run_id)

pal <- distinctColorPalette(n_distinct(long_filt$genus))

bind_rows(morsels) %>%
  ggplot(aes(x = genus, y = as.integer(n_samples), fill = genus)) +
  geom_boxplot(outlier.colour = NA) +
  scale_fill_manual(values = pal) +
  scale_y_continuous(breaks = seq(20)) +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3)) +
  labs(x = "Dominant genus", y = "No. patients")

bind_rows(morsels) %>%
  group_by(genus) %>%
  summarise(n = median(n_samples)) %>%
  arrange(desc(n))
long_filt %>%
  distinct(genus) %>% View()

long_filt %>%
  filter(genus == "Epilithonimonas")

species_df <- fread("results/metagenomic_out/RA.S.zeroed.csv")
species_df %>%
  pivot_longer(!run_id, names_to = "species", values_to = "rel_a") %>%
  filter(run_id == "46_2") %>%
  filter(rel_a > 0) %>%
  arrange(desc(rel_a))

long_df %>% filter(hap_vap_cap == "Water control") %>%
  filter(rel_a != 0) %>% 
  distinct(genus) %>% View()
