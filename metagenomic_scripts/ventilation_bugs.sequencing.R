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

df_filt <- fread("results/tax_classification_out/abundance_matrices/RA.G.zeroed.decontam.2.csv") %>%
  filter(run_id %in% meta_filt$run_id)

long_df <- df_filt %>%
  pivot_longer(!run_id, names_to = "taxa", values_to = "rel_a") %>%
  filter(rel_a != 0) %>%
  left_join(meta_filt %>% select(run_id, ventilation)) 

# long_df %>%
#   group_by(hap_vap_cap, taxa) %>%
#   summarise(n = n)
microbe_counts <- long_df %>%
  group_by(taxa) %>%
  summarise(n = n_distinct(run_id)) %>%
  arrange(desc(n))

microbe_subcounts <- long_df %>%
  group_by(taxa, ventilation) %>%
  summarise(n = n_distinct(run_id)) %>%
  pivot_wider(id_cols = taxa, names_from = ventilation, values_from = n) %>%
  mutate(n_HAP = replace_na(HAP, 0)) %>%
  mutate(n_VAP = replace_na(VAP, 0)) %>%
  select(-HAP, -VAP)

microbe_count_filt <- microbe_counts %>%
  filter(n > 10) %>%
  left_join(microbe_subcounts)

hap_counts <- long_df %>%
  group_by(ventilation) %>%
  summarise(n_total = n_distinct(run_id))

plot_df <- long_df %>%
  filter(taxa %in% microbe_count_filt$taxa) %>%
  group_by(taxa, ventilation) %>%
  summarise(n = n_distinct(run_id)) %>%
  left_join(hap_counts) %>%
  mutate(perc = n / n_total * 100)

# Plot differences
plot_df %>%
  mutate(taxa = factor(taxa, rev(microbe_counts$taxa))) %>%
  ggplot(aes(x = ventilation, y = taxa, fill = perc)) +
  geom_tile(color = "black") +
  scale_fill_viridis(option = "mako") +
  theme_classic() +
  geom_text(aes(label = signif(perc, 2)),
            color = "white") +
  labs(x = "Patient type", y = "Pathogen",
       fill = "% positive tests")

ggsave("results/metagenomic_out/ventilation_pathogens.sequencing.genus.pdf", dpi = 600,
       height = 6, width = 5)

n_mat <- plot_df %>%
  filter(taxa %in% microbe_count_filt$taxa) %>%
  pivot_wider(id_cols = taxa, names_from = ventilation, values_from = n) %>%
  mutate(HAP = replace_na(HAP, 0),
         VAP = replace_na(VAP, 0)) %>%
  column_to_rownames("taxa")

fisher.test(n_mat, simulate.p.value = T)

# Permutation test
observed <- plot_df %>%
  pivot_wider(id_cols = taxa, names_from = ventilation, values_from = perc) %>%
  mutate(obs_HAP = replace_na(HAP, 0),
         obs_VAP = replace_na(VAP, 0)) %>%
  select(-HAP, -VAP) %>%
  mutate(obs_ratio = obs_HAP / obs_VAP) %>%
  dplyr::select(taxa, obs_ratio) %>%
  left_join(microbe_count_filt) %>%
  mutate(parsed_taxa = str_glue("{taxa} (n={n})"))

set.seed(67)

n_iters <- 1000

perm_morsels <- foreach(i = seq(n_iters)) %do% {
  perm_labels <- meta_filt %>%
    mutate(ventilation = sample(meta_filt$ventilation, 
                                length(meta_filt$ventilation), 
                                replace = F)) %>%
    dplyr::select(run_id, ventilation)
  
  temp <- long_df %>%
    select(-ventilation) %>%
    left_join(perm_labels)
  
  temp_counts <- temp %>%
    group_by(ventilation) %>%
    summarise(n_total = n_distinct(run_id))
  
  temp_plot_df <- temp %>%
    filter(taxa %in% microbe_count_filt$taxa) %>%
    group_by(taxa, ventilation) %>%
    summarise(n = n_distinct(run_id)) %>%
    left_join(temp_counts) %>%
    mutate(perc = n / n_total) %>%
    pivot_wider(id_cols = taxa, names_from = ventilation, values_from = perc) %>%
    mutate(HAP = replace_na(HAP, 0),
           VAP = replace_na(VAP, 0)) %>%
    mutate(index = i)
  
  return(temp_plot_df)
}

perm_df <- bind_rows(perm_morsels) %>%
  mutate(ratio = HAP / VAP) %>%
  left_join(observed)

# Calculate p values
p_morsels <- foreach(bug_name = unique(perm_df$taxa)) %do% {
  temp <- perm_df %>%
    filter(taxa == bug_name)
  
  if(unique(temp$obs_ratio) > 1) {
    return(temp %>% 
             summarise(p = sum(ratio >= obs_ratio) / n_iters) %>%
             mutate(taxa = bug_name))
  } else {
    return(temp %>% 
             summarise(p = sum(ratio <= obs_ratio) / n_iters) %>%
             mutate(taxa = bug_name))
  }
}

# Filter tests
bug_filt <- observed %>%
  filter(n >= 30)
  # filter(obs_ratio > 1.5 | obs_ratio < 2/3)

pval_df <- bind_rows(p_morsels) %>%
  # filter(taxa %in% bug_filt$taxa) %>%
  mutate(adj_p = signif(p.adjust(p, "BH"), 3)) %>%
  arrange(adj_p) %>%
  mutate(significant = adj_p < 0.05) %>%
  left_join(observed) %>%
  mutate(parsed_taxa = str_glue("{taxa} (n={n})"))

pval_df %>%
  fwrite("results/metagenomic_out/ventilation_bug_permutation.genus.csv")

pval_df %>% 
  filter(significant) %>% View()
perm_df %>%
  left_join(microbe_count_filt %>% select(taxa, n)) %>%
  mutate(parsed_taxa = str_glue("{taxa} (n={n})")) %>%
  ggplot() +
  geom_density(aes(x = ratio)) +
  facet_wrap(.~parsed_taxa) +
  geom_vline(aes(xintercept = obs_ratio,
                 color = significant),
             lty = "dashed",
             data = pval_df) +
  geom_text(aes(x = obs_ratio, y = 1.5, 
                label = str_glue("adj. p={adj_p}"), 
                color = significant),
            data = pval_df,
            hjust = 0) +
  geom_text(aes(x = Inf, y = 1, label = str_glue("n={n}")),
            color = "black",
            data = microbe_count_filt,
            hjust = 0) +
  theme_bw() +
  labs(x = "Non-ventilated/ventilated", y = "Permutation density", color = "Adj. p<0.05")

ggsave("results/metagenomic_out/ventilation_bug_permutation.genus.pdf", dpi = 600,
       height = 8, width = 11)
