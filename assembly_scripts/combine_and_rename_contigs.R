rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(doParallel)
require(Biostrings)

file_dir <- "results/assembly_out/contig_out/"

file_list <- list.files(file_dir, full.names = T)

all <- foreach(file_name = file_list, .combine = "c") %do% {
  id <- gsub(file_dir, "", file_name)
  id <- gsub(".fna", "", id)
  fna <- readDNAStringSet(file_name)
  names(fna) <- str_glue("{id}.{names(fna)}")
  return(fna)
}

writeXStringSet(all, "results/assembly_out/all_contigs.renamed.fna")
