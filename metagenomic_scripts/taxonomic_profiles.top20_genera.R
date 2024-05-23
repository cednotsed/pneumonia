setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(vegan)
require(hillR)
require(ggpubr)
require(randomcoloR)

meta_filt <- fread("data/metadata/parsed_patient_metadata.filt.csv")
high_microbe <- fread("results/qc_out/high_microbe_samples.G.csv")

df_filt <- fread("results/metagenomic_out/RA.G.zeroed.csv") %>%
  filter(run_id %in% c(meta_filt$run_id))

long_df <- df_filt %>%
  pivot_longer(!run_id, names_to = "genus", values_to = "rel_a") %>%
  left_join(meta_filt)

meta_filt %>%
  filter(hap_vap_cap == "VAP") %>%
  summarise(n = n())

long_df %>%
  filter(hap_vap_cap == "VAP") %>%
  filter(rel_a > 0.5) %>%
  group_by(genus) %>%
  summarise(n_samples = n_distinct(run_id)) %>%
  arrange(desc(n_samples))

total_RA <- long_df %>%
  group_by(genus) %>%
  summarise(total_rel_a = sum(rel_a)) %>%
  arrange(desc(total_rel_a)) %>%
  head(30)

# Sample counts
count_df <- long_df %>%
  group_by(hap_vap_cap) %>%
  summarise(n = n_distinct(run_id))

# HCLUST
df_mat <- df_filt %>% column_to_rownames("run_id")
df_dist <- vegdist(df_mat, method = "bray")
clustering <- hclust(df_dist)

pal <- distinctColorPalette(n_distinct(long_df$genus))

long_df %>%
  filter(rel_a > 0.1) %>%
  mutate(run_id = factor(run_id, rownames(df_mat)[clustering$order])) %>%
  filter(genus %in% total_RA$genus) %>%
  ggplot(aes(x = run_id, y = rel_a, fill = genus)) +
  geom_bar(stat = "identity",
           position = "stack",
           color = "black") +
  scale_fill_manual(values = pal) +
  theme_bw() +
  facet_grid(cols = vars(hap_vap_cap), 
             space = "free", scale = "free") +
  theme(axis.text.x = element_blank(),
        # axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position = "top") +
  labs(x = "Patient", y = "Relative abundance",
       fill = "Culture/PCR-positive species") +
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
