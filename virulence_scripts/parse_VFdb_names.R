rm(list = ls())
setwd("C:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)

fna <- readDNAStringSet("databases/VFdb_170524/VFDB_setA_nt.170524.fna")
names(fna) <- str_split(names(fna), "\\(", simplify = T)[, 1]

writeXStringSet(fna, "databases/VFdb_170524/VFDB_setA_nt.170524.parsed_names.fna")
