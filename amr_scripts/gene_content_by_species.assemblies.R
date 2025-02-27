rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(ggpubr)
require(randomcoloR)
require(Hmisc)

meta <- fread("data/metadata/parsed_patient_metadata.filt.csv") %>%
  filter(high_microbe_count) %>%
  filter(hap_vap_cap %in% c("HAP", "VAP"))

# Get number of high quality contigs
checkm2 <- fread("results/assembly_out/assembly_merged_metadata.csv") %>%
  filter(run_id %in% meta$run_id) %>%
  filter(gtdbtk_warnings == "N/A") %>%
  filter(species != "") %>%
  filter(checkm2_contamination <= 5)

checkm2 %>%
  summarise(min = min(genome_size),
            max = max(genome_size))

# AMR content
merged <- fread("results/amr_out/amr_matrices/resfinder.assemblies.gene_content.csv") %>%
  filter(run_id %in% meta$run_id)

# Get class meta
class_meta <- merged %>%
  distinct(gene, class) %>%
  arrange(desc(gene), desc(class))

class_parsed <- foreach(gene_name = unique(class_meta$gene), .combine = "bind_rows") %do% {
  temp <- class_meta %>%
    filter(gene == gene_name)
  tibble(gene = gene_name, classes = paste0(temp$class, collapse = "\n"))
}

merged_filt <- merged %>%
  filter(run_id %in% meta$run_id) %>%
  filter(gtdbtk_warnings == "N/A") %>%
  filter(checkm2_contamination <= 5) %>%
  # Each run can have multiple bins matching same species
  distinct(run_id, species, gene, parsed_gene) %>%
  left_join(class_parsed) %>%
  mutate(gene_annot = str_glue("{parsed_gene}\n({classes})"))

sp_counts <- checkm2 %>%
  group_by(species) %>%
  summarise(n = n_distinct(run_id)) %>%
  arrange(desc(n))

sp_filt <- sp_counts %>%
  filter(species != "") %>%
  head(5)

plot_df <- merged_filt %>%
  filter(species %in% sp_filt$species) %>%
  mutate(species = factor(species, rev(sp_filt$species))) %>%
  group_by(species, gene_annot, parsed_gene) %>%
  summarise(n = n_distinct(run_id)) %>%
  arrange(desc(n)) 

plot_df %>%
  mutate(gene_annot = factor(gene_annot, unique(plot_df$gene_annot))) %>%
  ggplot(aes(x = gene_annot, y = species, fill = n)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low = "blue", high = "red") +
  geom_text(aes(label = n),
            color = "white") +
  labs(x = "Resistance gene", y = "", fill = "No. assemblies") +
  theme_classic() +
  theme(axis.text.y = element_text(face = "italic"))

ggsave("results/amr_out/gene_frequency_by_species.assemblies.with_annot.pdf", dpi = 600, width = 30, height = 3)

plot_df %>%
  mutate(parsed_gene = factor(parsed_gene, unique(plot_df$parsed_gene))) %>%
  ggplot(aes(x = parsed_gene, y = species, fill = n)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low = "blue", high = "red", limits = c(1, 12), breaks = c(3, 6, 9, 12)) +
  geom_text(aes(label = n),
            color = "white") +
  labs(x = "Resistance gene", y = "", fill = "No. assemblies") +
  theme_classic() +
  theme(axis.text.y = element_text(face = "italic"),
        axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1, face = "italic"))

ggsave("results/amr_out/gene_frequency_by_species.assemblies.pdf", dpi = 600, width = 12, height = 3)


# S. epidermidis
plot_df %>%
  filter(species == "Staphylococcus epidermidis") %>%
  ungroup() %>%
  distinct(parsed_gene)

base_count <- checkm2 %>%
  group_by(species) %>%
  summarise(sum_bases = sum(genome_size))

plot_df2 <- merged_filt %>%
  filter(species %in% sp_filt$species) %>%
  mutate(species = factor(species, rev(sp_filt$species))) %>%
  group_by(species, gene_annot, parsed_gene) %>%
  summarise(n = n_distinct(run_id)) %>% View()
  group_by(species) %>%
  summarise(n_types = n_distinct(parsed_gene),
            gene_freq = sum(n)) %>%
  ungroup() %>%
  left_join(base_count) %>%
  mutate(diversity_ratio = n_types / log10(sum_bases),
         freq_ratio = gene_freq / log10(sum_bases)) %>%
  arrange(desc(freq_ratio)) %>%
  pivot_longer(!c(species, n_types, gene_freq, sum_bases), names_to = "type", values_to = "ratio")

plot_df2 %>%
  mutate(species = factor(species, unique(plot_df2$species))) %>%
  ggplot(aes(x = species, y = ratio, fill = type)) +
  geom_col(color = "black") +
  facet_grid(rows = vars(type)) +
  theme_bw() +
  theme(legend.position = "none")

ggsave("results/amr_out/ratio_barplot.pdf", dpi = 600, width = 9, height = 4)
checkm2 %>%
  filter(checkm2_completeness >= 90) %>%
  summarise(median_n50 = median(contig_n50),
            n_genera = n_distinct(genus),
            n_species = n_distinct(species))
  nrow()
  

