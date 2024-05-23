rm(list = ls())
setwd("C:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)

path_list <- list.files("data/genomes/bug_references/", full.names = T)
file_list <- list.files("data/genomes/bug_references/", full.names = F)

df <- fread("data/metadata/bugs_to_download.curated.csv")
parsed <- tibble(file_path = paste0("/mnt/c/git_repos/pneumonia/", path_list),
       file_name = file_list) %>% 
  separate(file_name, c("acc1", "acc2"), "\\_") %>%
  mutate(accession = str_glue("{acc1}_{acc2}")) %>%
  left_join(df) %>%
  mutate(abundance = 1 / nrow(df) * 100) %>%
  select(bugs, file_path, abundance)

parsed %>%
  select(bugs, file_path) %>%
  fwrite("data/metadata/bugs.nanosim.meta.tsv",
         sep = "\t",
         col.names = F,
         eol = "\n")

parsed %>%
  select(Size = bugs, `1600` = abundance) %>%
  fwrite("data/metadata/bugs.nanosim.abundance.tsv",
         sep = "\t",
         eol = "\n")
read_df <- fread("results/metagenomic_out/read_counts.S.zeroed.csv")

read_df %>%
  pivot_longer(!run_id, names_to = "taxon", values_to = "read_count") %>%
  arrange(desc(read_count)) %>% View()


fread
