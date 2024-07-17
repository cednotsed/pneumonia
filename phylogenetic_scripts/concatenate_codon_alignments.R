rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)

bug_prefix <- "s_aureus"
file_list <- list.files(str_glue("results/pangenome_out/{bug_prefix}.filt/aligned_protein_sequences"), full.names = T)

assemblies <- deframe(fread(str_glue("data/metadata/bug_metadata/{bug_prefix}.assemblies.txt"),
                      header = F) %>%
  mutate(genome_name = gsub(".fna", "", V1)) %>%
  select(genome_name))

# Get intersecting genomes while keeping all assemblies
morsels <- foreach(file_name = file_list) %do% {
  # file_name <- file_list[1]
  temp <- readAAStringSet(file_name)
  parsed_names <- str_split(names(temp), ";", simplify = T)[, 1]
  
  # Number of assemblies carrying gene
  n_ass <- sum(assemblies %in% parsed_names)
  
  # Number of assemblies with duplicated genes
  n_dup_ass <- sum(duplicated(parsed_names[parsed_names %in% assemblies]))
  
  # Total number of genomes with duplicated genes
  n_dup <- sum(duplicated(parsed_names))
  
  tibble(file_name = file_name,
         n_assemblies = n_ass,
         n_assembly_duplicates = n_dup_ass,
         n_duplicates = n_dup)
}

merged_filt <- bind_rows(morsels) %>%
  filter(n_assemblies == length(assemblies)) %>%
  filter(n_assembly_duplicates == 0) %>%
  filter(n_duplicates == 0)

# Second pass to concatenate all alignments
morsels <- foreach(file_name = file_list) %do% {
  temp <- readAAStringSet(file_name)
  parsed_names <- str_split(names(temp), ";", simplify = T)[, 1]
  
  # Number of assemblies carrying gene
  n_ass <- sum(assemblies %in% parsed_names)
  
  # Number of assemblies with duplicated genes
  n_dup_ass <- sum(duplicated(parsed_names[parsed_names %in% assemblies]))
  
  # Total number of genomes with duplicated genes
  n_dup <- sum(duplicated(parsed_names))
  
  tibble(file_name = file_name,
         n_assemblies = n_ass,
         n_assembly_duplicates = n_dup_ass,
         n_duplicates = n_dup)
}