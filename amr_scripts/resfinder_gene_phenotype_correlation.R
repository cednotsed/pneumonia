rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(ggpubr)
require(randomcoloR)
require(igraph)

mat <- fread("results/amr_out/amr_matrices/resfinder.amr.decontam.csv") %>%
  column_to_rownames("run_id")

cor_mat <- cor(mat)

edgelist <- as.data.frame(cor_mat) %>%
  rownames_to_column("abx1") %>%
  pivot_longer(!abx1, names_to = "abx2", values_to = "corr") %>%
  filter(abx1 != abx2) %>%
  arrange(desc(corr)) %>%
  filter(corr > 0.90)
  
edgelist %>% 
  View()
edgelist %>%
  fwrite("results/amr_out/amr_matrices/resfinder.amr.decontam.gt90.edgelist.csv")
