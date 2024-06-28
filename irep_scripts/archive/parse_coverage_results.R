rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)

map_dir <- "results/irep_out/mapping_out/"

map_list <- list.files(map_dir, full.names = T)
map_list <- map_list[grepl("cov_stats", map_list)]

map_df <- foreach(file_name = map_list, .combine = "bind_rows") %do% {
  temp <- fread(file_name) %>%
    mutate(run = file_name)
  
  temp
}

map_parsed <- map_df %>%
  rename(chrom = `#rname`) %>%
  mutate(run = gsub(map_dir, "", run)) %>%
  mutate(run = gsub("INHALE_FRESH_|barcode|taxid", "", run)) %>%
  mutate(run = gsub("-0", "_", run)) %>%
  mutate(run = gsub("-", "_", run)) %>%
  separate(run, into = c("run_id", "taxid"), "\\.") %>%
  select(run_id, taxid, chrom, 
         genome_length = endpos, n_reads = numreads, coverage_breadth = coverage,
         mean_depth = meandepth)

map_parsed %>%
  fwrite("results/irep_out/coverage_results.parsed.tsv")

