rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)

irep_dir <- "results/irep_out/irep_results_out/"

irep_list <- list.files(irep_dir, full.names = T)
irep_list <- irep_list[grepl(".tsv", irep_list)]

irep_df <- foreach(file_name = irep_list, .combine = "bind_rows") %do% {
  temp <- fread(file_name)
  library_name <- colnames(temp)[4]
  
  temp_parsed <- temp %>%
    mutate(run = library_name) %>%
    mutate(across(everything(), as.character))
  
  colnames(temp_parsed) <- c("ref", "Ori", "Ter", "bPTR", "run")
  
  temp_parsed
}

irep_parsed <- irep_df %>%
  separate(ref, into = c(rep(NA, 9), "ref"), "/") %>% 
  separate(ref, into = c("acc1", "acc2"), "_") %>%
  mutate(ref = str_glue("{acc1}_{acc2}")) %>% 
  mutate(run = gsub("/SAN/ugi/HAP_VAP/pneumonia/results/irep_out/mapping_out/", "", run)) %>%
  mutate(run = gsub("INHALE_FRESH_|barcode|taxid", "", run)) %>%
  mutate(run = gsub("-0", "_", run)) %>%
  mutate(run = gsub("-", "_", run)) %>%
  separate(run, into = c("run_id", "taxid"), "\\.") %>%
  mutate(bPTR = as.numeric(bPTR)) %>%
  select(run_id, ref, taxid, Ori, Ter, bPTR)

fwrite(irep_parsed, "results/irep_out/irep_results.parsed.tsv")
