rm(list = ls())
setwd("C:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(ggpubr)
require(viridis)

meta <- fread("data/metadata/parsed_patient_metadata.filt.csv")

meta_filt <- meta %>%
  filter(hap_vap_cap %in% c("HAP", "VAP")) %>%
  filter(high_microbe_count) %>%
  mutate(ventilation = ifelse(ventilation, "VAP", "HAP"))

micro_meta <- fread("data/metadata/parsed_microbiology_results.csv") %>%
  filter(run_id %in% meta_filt$run_id) %>%
  filter(validity != "Invalid") %>%
  left_join(meta_filt %>% dplyr::select(run_id, hap_vap_cap, hap_vap2, ventilation)) %>%
  mutate(bugs = ifelse(grepl("Enterobacter", bugs), "Enterobacter spp.", bugs)) %>%
  mutate(bugs = ifelse(grepl("Acinetobacter", bugs), "Acinetobacter spp.", bugs)) %>%
  mutate(bugs = ifelse(grepl("Proteus", bugs), "Proteus spp.", bugs)) %>%
  filter(bugs != "Coliform") 

microbe_counts <- micro_meta %>%
  group_by(bugs) %>%
  summarise(n = n_distinct(run_id)) %>%
  arrange(desc(n))

microbe_count_filt <- microbe_counts %>%
  filter(n > 10) %>%
  filter(bugs != "Negative")

hap_counts <- micro_meta %>%
  group_by(ventilation) %>%
  summarise(n_total = n_distinct(run_id))

plot_df <- micro_meta %>%
  filter(bugs %in% microbe_count_filt$bugs) %>%
  group_by(bugs, ventilation) %>%
  summarise(n = n_distinct(run_id)) %>%
  left_join(hap_counts) %>%
  mutate(perc = n / n_total * 100)

# Permutation test
observed <- plot_df %>%
  pivot_wider(id_cols = bugs, names_from = ventilation, values_from = perc) %>%
  mutate(HAP = replace_na(HAP, 0),
         VAP = replace_na(VAP, 0)) %>%
  mutate(obs_ratio = HAP / VAP) %>%
  dplyr::select(bugs, obs_ratio)

set.seed(67)

n_iters <- 1000

perm_morsels <- foreach(i = seq(n_iters)) %do% {
  perm_labels <- meta_filt %>%
    mutate(ventilation = sample(meta_filt$ventilation, 
                                length(meta_filt$ventilation), 
                                replace = F)) %>%
    dplyr::select(run_id, ventilation)
  
  temp <- micro_meta %>%
    select(-ventilation) %>%
    left_join(perm_labels)
  
  temp_counts <- temp %>%
    group_by(ventilation) %>%
    summarise(n_total = n_distinct(run_id))
  
  temp_plot_df <- temp %>%
    filter(bugs %in% microbe_count_filt$bugs) %>%
    group_by(bugs, ventilation) %>%
    summarise(n = n_distinct(run_id)) %>%
    left_join(temp_counts) %>%
    mutate(perc = n / n_total) %>%
    pivot_wider(id_cols = bugs, names_from = ventilation, values_from = perc) %>%
    mutate(HAP = replace_na(HAP, 0),
           VAP = replace_na(VAP, 0)) %>%
    mutate(index = i)
  
  return(temp_plot_df)
}

perm_df <- bind_rows(perm_morsels) %>%
  mutate(ratio = HAP / VAP) %>%
  left_join(observed)

# Calculate p values
p_morsels <- foreach(bug_name = unique(perm_df$bugs)) %do% {
  temp <- perm_df %>%
    filter(bugs == bug_name)
  
  if(unique(temp$obs_ratio) > 1) {
    return(temp %>% 
             summarise(p = sum(ratio >= obs_ratio) / n_iters) %>%
             mutate(bugs = bug_name))
  } else {
    return(temp %>% 
             summarise(p = sum(ratio <= obs_ratio) / n_iters) %>%
             mutate(bugs = bug_name))
  }
}

# Filter tests
# bug_filt <- observed

pval_df <- bind_rows(p_morsels) %>%
  # filter(bugs %in% bug_filt$bugs) %>%
  mutate(adj_p = signif(p.adjust(p, "BH"), 3)) %>%
  mutate(significant = adj_p < 0.05) %>%
  arrange(adj_p)

pval_df %>%
  fwrite("results/metagenomic_out/ventilation_bug_permutation.PCR.csv")

perm_df %>%
  ggplot(aes(x = ratio)) +
  geom_density() +
  geom_vline(aes(xintercept = obs_ratio),
             color = "red",
             lty = "dashed",
             data = observed) +
  facet_wrap(.~bugs) +
  geom_text(aes(x = -Inf, y = 2, 
                label = str_glue("adj. p={adj_p}"), 
                color = significant),
            data = pval_df,
            hjust = 0) +
  geom_text(aes(x = Inf, y = 3, label = str_glue("n={n}")),
            color = "black",
            data = microbe_count_filt,
            hjust = 0)

# Plot differences
plot_df %>%
  mutate(pathogen_type = grepl("virus", bugs)) %>%
  mutate(bugs = factor(bugs, rev(microbe_counts$bugs))) %>%
  ggplot(aes(x = ventilation, y = bugs, fill = perc)) +
  geom_tile(color = "black") +
  scale_fill_viridis(option = "mako") +
  theme_classic() +
  geom_text(aes(label = signif(perc, 2)),
            color = "white") +
  facet_grid(rows = vars(pathogen_type), scales = "free", space = "free") +
  labs(x = "Patient type", y = "Pathogen",
       fill = "% positive tests")

ggsave("results/metagenomic_out/ventilation_pathogens.PCR.pdf", dpi = 600,
       height = 3, width = 5)

n_mat <- plot_df %>%
  filter(bugs %in% microbe_count_filt$bugs) %>%
  pivot_wider(id_cols = bugs, names_from = hap_vap_cap, values_from = n) %>%
  mutate(HAP = replace_na(HAP, 0),
         VAP = replace_na(VAP, 0)) %>%
  column_to_rownames("bugs")

fisher.test(n_mat, simulate.p.value = T)

