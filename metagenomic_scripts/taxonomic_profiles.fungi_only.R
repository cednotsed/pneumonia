setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)

meta <- fread("data/metadata/parsed_patient_metadata.csv")
taxa_meta <- fread("data/metadata/k2_database_taxonomy.csv")

fungi <- taxa_meta %>%
  filter(group == "Fungi")
rank <- "S"

RA_df <- fread(str_glue("results/metagenomic_out/RA.{rank}.zeroed.csv")) %>%
  # left_join(meta %>% select(run_id, hap_vap_cap)) %>%
  # filter(hap_vap_cap != "Water control") %>%
  # select(-hap_vap_cap) %>%
  column_to_rownames("run_id")

RA_filt <- RA_df %>%
  select(any_of(fungi$species))

long_df <- RA_filt %>%
  rownames_to_column("run_id") %>%
  pivot_longer(!run_id, names_to = "species", values_to = "rel_a") %>%
  left_join(meta)

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
ggsave("results/metagenomic_out/relative_abundance.fungi_only.pdf", width = 12, height = 5)


plot_df <- long_df %>%
  filter(rel_a > 0) %>%
  group_by(species) %>%
  summarise(n = n_distinct(run_id)) %>%
  arrange(desc(n))

plot_df %>%
  mutate(species = factor(species, unique(plot_df$species))) %>%
  ggplot(aes(x = species, y = n)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, 
                                   hjust = 1, 
                                   vjust = 0.3)) +
  labs(x = "Fungal species", y = "No. patients")

# SEPARATE GRAPHS
plt_list <- foreach(category = unique(long_df$hap_vap_cap)) %do% {
  long_df %>%
    filter(hap_vap_cap == category) %>%
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
              data = count_df %>% filter(hap_vap_cap == category),
              vjust = 1.1,
              hjust = -0.1)
}

pdf("results/metagenomic_out/relative_abundance.fungi_only.separate.pdf",
    width = 15, height = 5)
for(plt in plt_list) {
  print(plt)
}
dev.off()