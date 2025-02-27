rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(vegan)
require(ape)
require(gghalves)
require(ggnewscale)
require(paletteer)

file_list <- list.files("results/ML_out/results_out/", full.names = T)

morsels <- foreach(file_name = file_list) %do% {
  fread(file_name)
}

obs_df <- bind_rows(morsels) %>%
  select(g1, g2, threshold, test_precision, test_recall, test_F1, test_AUROC) %>%
  mutate(groups = str_glue("{g1} vs. {g2}")) %>%
  mutate(threshold = factor(threshold, c(1, 2, 5, 10))) %>%
  filter(g1 != "HAP")

perm_dir <- "results/ML_out/permutation_out/"
perm_list <- list.files(perm_dir, full.names = T)

perm_morsels <- foreach(perm_name = perm_list) %do% {
  file_suffix <- gsub(perm_dir, "", perm_name)
  suffix_vec <- str_split(file_suffix, "\\.")[[1]]
  g1 <- suffix_vec[1]
  g2 <- suffix_vec[2]
  threshold <- suffix_vec[3]
  
  fread(perm_name) %>%
    mutate(g1 = g1, g2 = g2, threshold = threshold)
}

perm_df <- bind_rows(perm_morsels) %>%
  mutate(groups = str_glue("{g1} vs. {g2}")) %>%
  mutate(threshold = factor(threshold, c(1, 2, 5, 10))) %>%
  filter(g1 != "HAP")

max_values <- perm_df %>%
  group_by(groups, threshold) %>%
  summarise(max_F1 = max(test_F1))

p_morsels <- foreach(i = seq(nrow(obs_df))) %do% {
  row <- obs_df[i, ]
  p_temp <- row %>%
    select(groups, threshold, obs_test_F1 = test_F1) %>%
    left_join(perm_df) %>%
    summarise(p = sum(test_F1 >= obs_test_F1) / 500)
  
  row %>%
    select(g1, g2, groups, threshold) %>%
    mutate(pval = p_temp$p)
}

p_df <- bind_rows(p_morsels) %>%
  mutate(adj_p = signif(p.adjust(pval, "BH"), 1)) %>%
  left_join(max_values) %>%
  mutate(threshold = factor(threshold, c(1, 2, 5, 10))) %>%
  filter(g1 != "HAP")

dodge_w <- 1
pal <- c("#E1E1E1FF", "#7D96AFFF", "#647D96FF", "#184860FF")

# Observed
obs_df %>%
  filter(threshold == 2) %>%
  select(groups, test_precision, test_recall, test_F1) %>%
  pivot_longer(!groups, names_to = "metric", values_to = "score") %>%
  ggplot(aes(x = groups, y = score, fill = metric)) +
  geom_col(position = position_dodge(width = dodge_w),
           color = "black") +
  scale_fill_manual(values = c("#E1E1E1FF", "#7D96AFFF", "#7F8C72FF")) +
  geom_text(aes(x = groups, y = max_F1 + 0.05, fill = threshold, label = str_glue("p={pval}")),
            position = position_dodge(width = dodge_w),
            angle = 90,
            data = p_df)


perm_df %>%
  ggplot() +
  geom_half_violin(aes(x = groups, y = test_F1, fill = threshold),
                   position = position_dodge(width = dodge_w)) +
  scale_fill_manual(values = pal) +
  new_scale_fill() +
  geom_point(aes(x = groups, y = test_F1, groups = threshold),
             fill = "#D5B77DFF",
             color = "black",
             position = position_dodge(width = dodge_w),
             size = 3,
             pch = 23,
             data = obs_df) +
  scale_fill_discrete(guide = "none") +
  geom_text(aes(x = groups, y = max_F1 + 0.05, fill = threshold, label = str_glue("p={pval}")),
            position = position_dodge(width = dodge_w),
            angle = 90,
            data = p_df) +
  theme_classic() +
  labs(x = "Comparison", y = "Test F1 (5-fold CV)")

ggsave("results/ML_out/permutation_test.pdf", width = 5, height = 4, dpi = 600)
