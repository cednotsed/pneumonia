setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(Biostrings)

fna <- readDNAStringSet("data/genomes/all_irep_refs.fna")
fna_filt <- fna[!grepl("plasmid", names(fna))]


writeXStringSet(fna_filt, "data/genomes/all_irep_refs.no_plasmids.fna")
## REMEMBER TO DOS2UNIX
