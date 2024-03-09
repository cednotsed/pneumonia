setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(vegan)
require(hillR)
require(ggpubr)
require(randomcoloR)

meta <- fread("data/metadata/parsed_patient_metadata.filt.csv")
high_microbe <- fread("results/qc_out/high_microbe_samples.S.csv")

df_filt <- fread("results/metagenomic_out/RA.S.zeroed.csv")

# Get bug list
bug_list <- unique(c(meta$micro_organism, meta$biofire_organism, meta$curetis_organism))
bug_list <- gsub("\\s*\\([^\\)]+\\)","", bug_list)
bug_list <- unique(unlist(str_split(bug_list, ";")))
bug_list <- bug_list[!(bug_list %in% c("Negative", "", "Invalid test"))]
# Get species only
bug_list <- bug_list[!grepl("complex|group|Coliform|spp|virus|-|Pneumocystis", bug_list)]
# bug_list <- bug_list[bug_list %in% c("Coliform", "Coronavirus", "Adenovirus")]

long_df <- df_filt %>%
  select(all_of(c("run_id", bug_list))) %>%
  pivot_longer(!run_id, names_to = "species", values_to = "rel_a") %>%
  left_join(meta)

long_df
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

pdf("results/metagenomic_out/relative_abundance.positive_bugs.separate.pdf",
    width = 15, height = 5)
  for(plt in plt_list) {
    print(plt)
  }
dev.off()

