rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)

tax_meta <- fread("databases/k2_pluspfp_20240112/inspect.txt")
fungi <- tax_meta[501:801, ]$V6

meta_filt <- fread("data/metadata/parsed_patient_metadata.filt.csv") %>%
  filter(high_microbe_count) %>%
  filter(hap_vap_cap %in% c("HAP", "VAP"))

df_filt <- fread("results/tax_classification_out/abundance_matrices/RA.G.zeroed.decontam.2.csv") %>%
  filter(run_id %in% meta_filt$run_id)

long_df <- df_filt %>%
  pivot_longer(!run_id, names_to = "taxa", values_to = "rel_a") %>%
  filter(rel_a != 0)

fungal_df <- long_df %>%
  filter(taxa %in% fungi)

# PCR and Culture of fungal patients
micro_meta <- fread("data/metadata/parsed_microbiology_results.csv") %>%
  filter(bugs != "Invalid test")

# Plot fungal abundance
# fungal_mat <- fungal_df %>%
#   pivot_wider(id_cols = run_id, names_from = taxa, values_from = rel_a) %>%
#   column_to_rownames("run_id")

plot_df <- fungal_df %>%
  filter(rel_a != 0) %>%
  left_join(meta_filt %>% select(run_id, hap_vap_cap))

fungi_filt <- plot_df %>% 
  group_by(taxa) %>%
  summarise(n = n()) %>%
  filter(n != 1)

order_df <- plot_df %>%
  mutate(taxa %in% fungi_filt$taxa) %>%
  group_by(run_id) %>%
  summarise(total_rel_a = sum(rel_a)) %>%
  arrange(desc(total_rel_a))

count_df <- plot_df %>% 
  filter(rel_a != 0) %>%
  filter(taxa %in% fungi_filt$taxa) %>%
  distinct(run_id, hap_vap_cap) %>%
  group_by(hap_vap_cap) %>%
  summarise(n = n())

plot_df %>% 
  filter(rel_a != 0) %>% 
  group_by(hap_vap_cap) %>% summarise(n = n_distinct(run_id))

# pal <- distinctColorPalette(n_distinct(fungi_filt$taxa))
pal <- c("#DB8F95", "#A4CBD0", "black", "#BA6AD9")

plot_df %>%
  filter(rel_a != 0) %>%
  filter(taxa %in% fungi_filt$taxa) %>%
  mutate(run_id = factor(run_id, order_df$run_id)) %>%
  ggplot(aes(x = run_id, y = rel_a, fill = taxa)) +
  geom_bar(stat = "identity",
           position = "stack",
           color = "black") +
  scale_fill_manual(values = pal, na.translate = F) +
  theme_classic() +
  facet_grid(cols = vars(hap_vap_cap), scales = "free", space = "free") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "top") +
  labs(x = "Patient", y = "Relative abundance",
       fill = "Fungal genus") +
geom_text(aes(x = 0, y = Inf, 
              fill = NA, 
              label = str_glue("n={n}")),
          data = count_df,
          vjust = 1.1,
          hjust = -0.1)

ggsave("results/metagenomic_out/fungal.pdf", dpi = 600, width = 6, height = 4)

# Relative abundance of fungi
plot_df %>%
  filter(rel_a > 0.5) %>%
  distinct(run_id) %>%
  nrow()

# Controls
control_meta <- fread("data/metadata/parsed_patient_metadata.filt.csv") %>%
  filter(high_microbe_count) %>%
  filter(hap_vap_cap == "healthy")

control_RA <- fread("results/tax_classification_out/abundance_matrices/RA.G.zeroed.csv") %>%
  filter(run_id %in% control_meta$run_id)

control_filt <- control_RA %>%
  pivot_longer(!run_id, names_to = "taxa", values_to = "rel_a") %>%
  filter(rel_a != 0)

control_filt %>% 
  filter(taxa %in% fungi)

# Undiagnosed
undiagnosed_ids <- fread("results/benchmarking_out/undiagnosed_ids.txt", header = F)$V1

plot_df %>%
  filter(run_id %in% undiagnosed_ids) %>%
  filter(rel_a > 0.5) %>%
  distinct(run_id) %>%
  nrow()

# Read
test_positivity <- fungal_df %>%
  distinct(run_id) %>%
  left_join(micro_meta) %>%
  group_by(run_id) %>%
  summarise(any_positive = any(bugs != "Negative"))

fungal_df %>%
  left_join(test_positivity) %>%
  ggplot(aes(x = any_positive, y = rel_a)) +
  geom_boxplot() +
  geom_pwc()
