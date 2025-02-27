rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(ggpubr)
require(randomcoloR)
require(Hmisc)

meta <- fread("data/metadata/parsed_patient_metadata.filt.csv")
merged <- fread("results/amr_out/amr_matrices/resfinder.gene_content.csv") %>%
  filter(run_id %in% meta$run_id)

class_count <- merged %>%
  group_by(gene) %>% 
  summarise(n = n_distinct(class)) %>%
  arrange(desc(n)) %>%
  filter(n == 1)

class_meta <- merged %>%
  distinct(gene, class) %>%
  mutate(class = ifelse(gene %in% class_count$gene, class, "Multi-class")) %>%
  distinct()

plot_df <- merged %>%
  group_by(gene) %>% 
  summarise(n = n_distinct(run_id)) %>%
  arrange(desc(n)) %>%
  left_join(class_meta) %>%
  head(20)

pal <- c("cyan4", "darkolivegreen4", "darkgoldenrod3", 
         "#2e4057", "#d1495b", "cornflowerblue", 
         "darkorchid", "aquamarine", "antiquewhite")
plot_df %>%
  mutate(gene = factor(gene, plot_df$gene)) %>%
  ggplot(aes(x = gene, y = n / 239, fill = class)) +
  geom_col(color = "black") +
  scale_fill_manual(values = pal) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(x = "Resistance gene", y = "Prop. patients")


# ggsave("results/amr_out/mlsb_piechart.pdf", width = 5, height = 4, dpi = 600)

test <- fread("results/amr_out/resfinder_out/resfinder.parsed.decontam.csv") %>%
  filter(run_id %in% meta$run_id)

test %>%
  filter(wgs_predicted_phenotype == "Resistant") %>%
  filter(grepl("aph", genetic_background)) %>%
  distinct(antimicrobial)
  distinct(run_id)
