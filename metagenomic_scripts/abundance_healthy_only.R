setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(vegan)
require(hillR)
require(ggpubr)
require(randomcoloR)

meta_filt <- fread("data/metadata/parsed_patient_metadata.filt.csv")

df_filt <- fread("results/tax_classification_out/abundance_matrices/RA.G.zeroed.decontam.2.csv")

long_df <- df_filt %>%
  pivot_longer(!run_id, names_to = "genus", values_to = "rel_a") %>%
  left_join(meta_filt)

healthy <- long_df %>% 
  filter(hap_vap_cap == "healthy") %>%
  filter(rel_a != 0) %>%
  distinct(genus) 

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
  filter(hap_vap_cap == "healthy") %>%
  mutate(run_id = factor(run_id, rownames(df_mat)[clustering$order])) %>%
  filter(genus %in% healthy$genus) %>%
  ggplot(aes(x = run_id, y = rel_a, fill = genus)) +
  geom_bar(stat = "identity",
           position = "stack",
           color = "black") +
  scale_fill_manual(values = pal) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        # axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position = "top") +
  labs(x = "Patient", y = "Relative abundance",
       fill = "Microbial genus") +
  geom_text(x = 0, y = 1, 
            label = str_glue("n={n_distinct(healthy$run_id)}"))

ggsave("results/metagenomic_out/relative_abundance_healthy.pdf", width = 12, height = 5)
