rm(list = ls())
setwd("C:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)
require(ggpubr)
require(randomcoloR)

# Metadata
meta <- fread("data/metadata/parsed_patient_metadata.filt.csv")

meta_filt <- meta %>%
  filter(hap_vap_cap %in% c("HAP", "VAP"))

undiagnosed_ids <- fread("results/benchmarking_out/undiagnosed_ids.txt")$run_id

control_meta <- meta %>%
  filter(hap_vap_cap == "Water control")

micro_meta <- fread("data/metadata/parsed_microbiology_results.csv") %>%
  filter(run_id %in% undiagnosed_ids) %>%
  filter(run_id %in% meta_filt$run_id)

bact_meta <- fread("data/metadata/parsed_microbiology_results.bacterial_sp_only.csv") %>%
  filter(species_resolved) %>%
  filter(run_id %in% undiagnosed_ids) %>%
  filter(run_id %in% meta_filt$run_id)

# Sequencing detection
RA_df <- fread("results/tax_classification_out/abundance_matrices/RA.S.zeroed.decontam.2.csv")
read_df <- fread("results/tax_classification_out/abundance_matrices/read_counts.S.zeroed.decontam.2.csv")

long_RA_df <- RA_df %>%
  pivot_longer(!run_id, names_to = "species", values_to = "rel_a")

long_read_df <- read_df %>%
  pivot_longer(!run_id, names_to = "species", values_to = "read_count")

merged_df <- long_RA_df %>%
  filter(run_id %in% undiagnosed_ids) %>% 
  left_join(long_read_df %>% select(run_id, species, read_count)) %>%
  filter(rel_a != 0)

# Get dominant pathogen
morsels <- foreach(run_name = unique(merged_df$run_id)) %do% {
  merged_df %>%
    filter(run_id == run_name) %>%
    arrange(desc(rel_a)) %>%
    head(1)
}

plot_df <- bind_rows(morsels) %>%
  arrange(species)

pal <- distinctColorPalette(n_distinct(plot_df$species))
pal <- setNames(pal, unique(plot_df$species))

p1 <- plot_df %>%
  mutate(run_id = factor(run_id, unique(plot_df$run_id)),
         species = factor(species, unique(plot_df$species))) %>%
  ggplot(aes(x = run_id, y = rel_a, fill = species)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = pal) +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_text(face = "italic", size = 8),
        legend.title = element_text(angle = 90)) +
  geom_text(aes(label = read_count), vjust = 0.5, hjust = 1.1, size = 3, angle = 90) +
  labs(x = "Patient", y = "Relative abundance", fill = "Dominant species")

p2 <- plot_df %>%
  mutate(run_id = factor(run_id, unique(plot_df$run_id)),
         species = factor(species, unique(plot_df$species))) %>%
  ggplot(aes(y = log10(read_count))) +
  geom_boxplot(fill = "grey") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(y = "Log10(read count)")

p3 <- plot_df %>%
  mutate(run_id = factor(run_id, unique(plot_df$run_id)),
         species = factor(species, unique(plot_df$species))) %>%
  ggplot(aes(y = rel_a)) +
  geom_boxplot(fill = "grey") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(y = "Relative abundance")

ggarrange(p1, p2, p3,
          nrow = 1, ncol = 3,
          align = "v", 
          widths = c(5, 1, 1)) 

ggsave("results/benchmarking_out/undiagnosed_dominant_bugs.pdf", dpi = 600,
       width = 12, height = 4)

plot_df %>%
  summarise(median_read = median(read_count),
            median_a = median(rel_a))
# Save species
# plot_df %>%
#   distinct(run_id, species) %>%
#   fwrite("results/benchmarking_out/undiagnosed_dominant_bugs.csv")

# long_read_df %>%
#   filter(run_id %in% undiagnosed_ids) %>%
#   filter(!grepl("Candida|Aspergillus|Penicillium|virus|Talaromyces|Nakaseomyces|Saccharomyces", species)) %>% 
#   filter(read_count >= 10000) %>%
#   fwrite("data/metadata/irep/undiagnosed_bugs.csv", 
#          eol = "\n")
  
# Check controls
contam_RA_df <- fread("results/tax_classification_out/abundance_matrices/RA.S.zeroed.csv") %>%
  pivot_longer(!run_id, names_to = "species", values_to = "rel_a")

control_RA <- contam_RA_df %>%
  # filter(run_id %in% control_meta$run_id) %>%
  right_join(control_meta %>% select(run_id, run))

control_RA %>%
  group_by(run) %>%
  summarise(n = n_distinct(run_id)) %>%
  arrange(desc(n))

plot_df %>%
  left_join(meta %>%
              select(run, run_id)) %>%
  left_join(control_RA %>% select(run, species, control_rel_a = rel_a)) %>%
  mutate(fold_diff = rel_a / control_rel_a) %>% View()

