rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(ggpubr)
require(viridis)
require(Hmisc)

meta <- fread("data/metadata/parsed_patient_metadata.filt.csv") %>%
  filter(high_microbe_count)

df <- fread("results/amr_out/resfinder_out/resfinder.parsed.decontam.csv") %>%
  filter(hap_vap_cap %in% c("HAP", "VAP")) %>%
  mutate(class = capitalize(class))

count_df <- df %>%
  group_by(hap_vap_cap) %>%
  summarise(n_total = n_distinct(run_id))

abx_counts <- df %>%
  filter(wgs_predicted_phenotype == "Resistant") %>%
  group_by(class) %>%
  summarise(n_abx = n_distinct(run_id)) %>%
  filter(n_abx >= 10)

plot_df <- df %>%
  filter(wgs_predicted_phenotype == "Resistant") %>%
  filter(class %in% abx_counts$class) %>%
  group_by(hap_vap_cap, class) %>%
  summarise(n = n_distinct(run_id)) %>% 
  left_join(count_df) %>%
  mutate(perc = n / n_total * 100)

plot_df %>%
  ggplot(aes(x = hap_vap_cap, y = class, fill = perc)) +
  geom_tile(color = "black") +
  geom_text(aes(label = signif(perc, 2)), color = "white") +
  scale_fill_viridis(option = "cividis") +
  theme_classic() +
  labs(x = "Patient type", y = "Antimicrobial class", fill = "% patients")

ggsave("results/amr_out/amr_profiles.pdf", dpi = 600, width = 5, height = 3)

# Permutation test
observed <- plot_df %>%
  pivot_wider(id_cols = class, names_from = hap_vap_cap, values_from = perc) %>%
  mutate(HAP = replace_na(HAP, 0),
         VAP = replace_na(VAP, 0)) %>%
  mutate(obs_ratio = HAP / (VAP)) %>%
  dplyr::select(class, obs_ratio)

n_iters <- 1000
set.seed(66)

perm_morsels <- foreach(i = seq(n_iters)) %do% {
  perm_labels <- meta %>% 
    filter(run_id %in% df$run_id) %>%
    select(run_id, hap_vap_cap)
  
  perm_labels <- perm_labels %>%
    mutate(hap_vap_cap = sample(perm_labels$hap_vap_cap, length(perm_labels$hap_vap_cap), replace = F))

  temp <- df %>%
    select(-hap_vap_cap) %>%
    filter(class %in% abx_counts$class) %>%
    left_join(perm_labels) 
    
  temp_counts <- temp %>%
    group_by(hap_vap_cap) %>%
    summarise(n_total = n_distinct(run_id))
  
  temp %>%
    filter(wgs_predicted_phenotype == "Resistant") %>%
    group_by(hap_vap_cap, class) %>%
    summarise(n = n_distinct(run_id)) %>% 
    left_join(temp_counts) %>%
    mutate(perc = n / n_total * 100) %>%
    pivot_wider(id_cols = class, names_from = hap_vap_cap, values_from = perc) %>%
    mutate(HAP = replace_na(HAP, 0),
           VAP = replace_na(VAP, 0)) %>%
    mutate(index = i)
}

perm_df <- bind_rows(perm_morsels) %>%
  mutate(ratio = HAP / VAP) %>%
  left_join(observed)

# Calculate p values
p_morsels <- foreach(abx = unique(perm_df$class)) %do% {
  temp <- perm_df %>%
    filter(class == abx)
  
  if(unique(temp$obs_ratio) > 1) {
    return(temp %>% 
             summarise(p = sum(ratio >= obs_ratio) / n_iters) %>%
             mutate(class = abx))
  } else {
    return(temp %>% 
             summarise(p = sum(ratio <= obs_ratio) / n_iters) %>%
             mutate(class = abx))
  }
}

pval_df <- bind_rows(p_morsels) %>%
  mutate(adj_p = signif(p.adjust(p, "BH"), 3)) %>%
  arrange(adj_p) %>%
  mutate(significant = adj_p < 0.05) 

perm_df %>%
  ggplot(aes(x = ratio)) +
  geom_density() +
  geom_vline(aes(xintercept = obs_ratio),
             color = "red",
             lty = "dashed",
             data = observed) +
  facet_wrap(.~class) +
  geom_text(aes(x = -Inf, y = 2, 
                label = str_glue("adj. p={adj_p}"), 
                color = significant),
            data = pval_df,
            hjust = 0) +
  geom_text(aes(x = -Inf, y = 3, label = str_glue("n={n_abx}")),
            color = "black",
            data = abx_counts,
            hjust = 0)
# parsed %>%
#   filter(wgs_predicted_phenotype == "Resistant") %>%
#   filter(hap_vap_cap == "Water control") %>%
#   ggplot(aes(x = run_id, y = antimicrobial, fill = wgs_predicted_phenotype)) +
#   geom_tile() +
#   scale_fill_manual(values = c("steelblue", "red")) 
# distinct(antimicrobial)
# drug_list <- deframe(parsed %>%
#                        distinct(antimicrobial) %>%
#                        arrange(desc(antimicrobial)))
# parsed %>%
#   mutate(antimicrobial = factor(antimicrobial, drug_list)) %>%
#   ggplot(aes(x = run_id, y = antimicrobial, fill = wgs_predicted_phenotype)) +
#   geom_tile() +
#   scale_fill_manual(values = c("steelblue", "red")) +
#   facet_grid(cols = vars(hap_vap_cap), scales = "free", space = "free") +
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank()) +
#   labs(x = "Sample", y = "Antimicrobial",
#        fill = "Predicted phenotype")
# 
# ggsave("results/amr_out/amr_profiles.pdf", width = 15, height = 15)
# 
# parsed %>%
#   group_by(run_id) %>%
#   filter(wgs_predicted_phenotype == "Resistant") %>%
#   summarise(n = n_distinct(class)) %>%
#   ggplot(aes(x = n)) +
#   geom_histogram() +
#   labs(x = "No. antimicrobial resistance classes", 
#        y = "No. patients")
# 
# plot_df <- parsed %>%
#   group_by(antimicrobial) %>%
#   summarise(prop = sum(wgs_predicted_phenotype == "Resistant") / n()) %>%
#   arrange(desc(prop)) %>% 
#   head(20) %>%
#   mutate(antimicrobial = Hmisc::capitalize(antimicrobial))
# pal <- distinctColorPalette(n_distinct(plot_df$antimicrobial))
# plot_df %>%
#   filter(antimicrobial %in% plot_df$antimicrobial) %>%
#   mutate(antimicrobial = factor(antimicrobial, plot_df$antimicrobial),
#          prop = round(prop, 2)) %>%
#   ggplot(aes(x = antimicrobial, y = prop, fill = antimicrobial)) +
#   geom_bar(stat = "identity",
#            color = "black") +
#   # geom_text(aes(label = prop)) +
#   scale_fill_manual(values = pal) + 
#   labs(x = "Antimicrobial", y = "Prop. of patients resistant to") +
#   theme_classic() + 
#   theme(legend.position = "none",
#         axis.text.x = element_text(angle = 45, hjust = 1))
# 
# ggsave("results/amr_out/amr_barplots.pdf", dpi = 600, width = 5, height = 5)  
# ggsave("results/amr_out/amr_barplots.png", dpi = 600, width = 5, height = 5)  
# filter(wgs_predicted_phenotype == "Resistant") %>%
#   distinct(run_id)
# 
# summarise(n = n_distinct(antimicrobial))
# id_parsed <- gsub("INHALE_FRESH_|.pheno_table.txt|barcode", "", id_list)
# id_parsed <- gsub("-0", "-", id_parsed)
# id_parsed <- gsub("-", "_", id_parsed)
# id_parsed <- gsub("12a", "12", id_parsed)
# 
# id_parsed
# left_join(meta %>% select(run_id, hap_vap_cap)) 
# meta %>%
#   filter(run_id %in% id_parsed)
# parsed %>%
#   distinct()
# 
