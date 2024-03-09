rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)

all <- fread("data/metadata/irep/irep_pairs.txt",
             header = F) %>%
  # select(V1, V5) %>%
  mutate(run = gsub(".fastq.gz", "", V1)) %>%
  mutate(string = paste0(run, V5))
all
done <- fread("results/irep_out/done.txt",
              header = F) %>%
  mutate(V2 = gsub("taxid", "", V2)) %>%
  mutate(string = paste0(V1, V2))

missing <- all$string[!(all$string %in% done$string)]

all %>%
  filter(string %in% missing) %>%
  select(-run, -string) %>%
  fwrite("data/metadata/irep/missing_irep_pairs.txt",
         col.names = F,
         eol = "\n")
  




