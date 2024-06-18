rm(list = ls())
setwd("/flask/scratch/matthewsp/pneumonia")
require(tidyverse)
require(data.table)
require(foreach)

irep_dir <- "results/irep_out/irep_out/"

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
    mutate(ref = gsub("/flask/scratch/matthewsp/pneumonia/data/genomes/irep_filt/", "", ref)) %>%
    mutate(run = gsub("/flask/scratch/matthewsp/pneumonia/results/irep_out/mapping_out/|.sam", "", run))

fwrite(irep_df, "results/irep_out/irep_results.tsv")
