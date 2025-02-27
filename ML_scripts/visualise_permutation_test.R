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
  filter(g1 != "HAP") %>%
  filter(threshold == 2)

perm_dir <- "results/ML_out/permutation_out/"
perm_list <- list.files(perm_dir, full.names = T)
perm_list <- perm_list[grepl("\\.2\\.", perm_list)]

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
  filter(g1 != "HAP")

p_morsels <- foreach(i = seq(nrow(obs_df))) %do% {
  row <- obs_df[i, ]
  p_temp <- row %>%
    select(groups, 
           obs_test_F1 = test_F1, 
           obs_test_precision = test_precision,
           obs_test_recall = test_recall) %>%
    left_join(perm_df) %>%
    summarise(p_F1 = sum(test_F1 >= obs_test_F1) / 1000,
              p_precision = sum(test_precision >= obs_test_precision) / 1000,
              p_recall = sum(test_recall >= obs_test_recall) / 1000)
  
  row %>%
    select(g1, g2, groups) %>%
    bind_cols(p_temp)
}

max_values <- perm_df %>%
  select(groups, test_precision, test_recall, test_F1) %>%
  pivot_longer(!groups, names_to = "metric", values_to = "score") %>%
  group_by(groups, metric) %>%
  summarise(max_score = max(score))

p_df <- bind_rows(p_morsels) %>%
  select(groups, p_F1, p_precision, p_recall) %>%
  pivot_longer(!groups, names_to = "metric", values_to = "pval") %>%
  mutate(metric = gsub("p_", "test_", metric)) %>%
  left_join(max_values)

dodge_w <- 1
pal <- c("#E1E1E1FF", "#7D96AFFF", "#647D96FF", "#184860FF")

# Plot
plot_df <- obs_df %>%
  select(groups, test_precision, test_recall, test_F1) %>%
  pivot_longer(!groups, names_to = "metric", values_to = "score")

perm_df %>%
  select(groups, test_recall, test_F1, test_precision) %>%
  pivot_longer(!groups, names_to = "metric", values_to = "score") %>%
  ggplot() +
  geom_violin(aes(x = groups, y = score, fill = metric),
              position = position_dodge(width = dodge_w)) +
  scale_fill_manual(values = c("#E1E1E1FF", "#7D96AFFF", "#7F8C72FF")) +
  geom_point(aes(x = groups, y = score, groups = metric),
             position = position_dodge(width = dodge_w),
             color = "black",
             fill = "black",
             size = 3,
             pch = 23,
             data = plot_df) +
  geom_text(aes(x = groups, 
                y = max_score + 0.2, 
                label = str_glue("p={pval}"),
                groups = metric),
            position = position_dodge(width = dodge_w),
            angle = 90,
            data = p_df) +
  ylim(0, 1) +
  theme_classic() +
  theme(legend.position = "none")
  # geom_col(position = position_dodge(width = dodge_w),
  #          color = "black") +
  # scale_fill_manual(values = c("#E1E1E1FF", "#7D96AFFF", "#7F8C72FF")) +
  # geom_text(aes(x = groups, y = max_score + 0.05, label = str_glue("p={pval}")),
  #           position = position_dodge(width = dodge_w),
  #           angle = 90,
  #           data = p_df) +

ggsave("results/ML_out/permutation_test.pdf", width = 5, height = 3, dpi = 600)

