setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(vegan)
require(hillR)
require(ggpubr)
require(randomcoloR)

meta <- fread("data/metadata/parsed_patient_metadata.csv")
high_microbe <- fread("results/qc_out/high_microbe_samples.S.csv")

df_filt <- fread("results/metagenomic_out/RA.S.zeroed.csv")

long_df <- df_filt %>%
  pivot_longer(!run_id, names_to = "species", values_to = "rel_a") %>%
  left_join(meta)

long_df %>%
  group_by(run_id, hap_vap_cap) %>%
  summarise(max_a = max(rel_a)) %>%
  ggplot(aes(x = max_a, fill = hap_vap_cap)) +
  geom_histogram() +
  facet_grid(rows = vars(hap_vap_cap)) +
  scale_fill_manual(values = c("indianred", "blue", "goldenrod", "grey50")) +
  labs(x = "Max. species abundance", y = "No. samples",
       fill = "Sample type")

total_RA <- long_df %>%
  group_by(run_id) %>%
  summarise(total_rel_a = sum(rel_a)) %>%
  arrange(desc(total_rel_a))

# Sample counts
count_df <- long_df %>%
  group_by(hap_vap_cap) %>%
  summarise(n = n_distinct(run_id))

pal <- distinctColorPalette(n_distinct(long_df$species))

long_df %>%
  mutate(run_id = factor(run_id, unique(total_RA$run_id))) %>%
  ggplot(aes(x = run_id, y = rel_a, fill = species)) +
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

