rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)

meta <- fread("data/metadata/irep/irep.config.tsv")
irep_dir <- "results/irep_out/irep_results/"

irep_list <- list.files(irep_dir, full.names = T)
irep_list <- irep_list[grepl(".tsv", irep_list)]

#irep_list <- irep_list[1:40000]

morsels <- foreach(file_name = irep_list) %do% {
#  file_name = irep_list[49450]
  temp <- fread(file_name)

  if (nrow(temp) > 0) {
      library_name <- colnames(temp)[4]
      print(library_name)
    
    temp_parsed <- temp %>%
        mutate(run = library_name) %>%
        mutate(across(everything(), as.character))
  
    colnames(temp_parsed) <- c("ref", "Ori", "Ter", "bPTR", "run")
    return(temp_parsed)
  } else {
    return(NULL)
  }
}

irep_df <- bind_rows(morsels) %>%
    mutate(run = gsub("/mnt/c/git_repos/pneumonia/results/irep_out/mapping_out/", "", run)) %>%
    separate(run, c("run_id", "taxid"), "\\.") %>%
    left_join(meta %>% 
                select(taxa = V2, taxid = V3) %>%
                mutate(taxid = as.character(taxid)) %>%
                distinct())

irep_df %>%
  filter(bPTR != "n/a") %>% View()
  distinct(taxa)

fwrite(irep_df, "results/irep_out/irep_results.tsv")
