rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)

all <- fread("data/metadata/irep/irep_pairs.txt",
             header = F)

multi <- fread("data/metadata/irep/multi_chromosome_taxa.txt")
all %>%
  filter(V5 %in% multi$V1) %>%
  fwrite("data/metadata/irep/multi_chromosome_pairs.txt",
         col.names = F,
         eol = "\n")





