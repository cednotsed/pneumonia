rm(list = ls())
setwd("C:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)

meta <- fread("data/metadata/parsed_patient_metadata.filt.csv") %>%
  filter(hap_vap_cap %in% c("HAP", "VAP"))

morsels <- foreach(id = unique(meta$run_id)) %do% {
  print(id)
  temp <- meta %>%
    filter(run_id == id)
  
  micro <- tibble(run_id = id, 
                  bugs = str_split(temp$culture_organism, "\\;", simplify = T)[1, ],
                  method = "culture")
  biofire <- tibble(run_id = id, 
                  bugs = str_split(temp$biofire_organism, "\\;", simplify = T)[1, ],
                  method = "biofire", 
                  validity = temp$biofire_valid)
  curetis <- tibble(run_id = id, 
                    bugs = str_split(temp$curetis_organism, "\\;", simplify = T)[1, ],
                    method = "curetis", 
                    validity = temp$curetis_valid)
  
  bind_rows(micro, biofire, curetis)
  
}

parsed <- bind_rows(morsels) %>% 
  filter(bugs != "") %>%
  mutate(bugs = gsub("\\(|>=10\\^7\\)|\\+|\\)|\\^", "", bugs)) %>%
  mutate(bugs = gsub("[[:digit:]]+", "", bugs)) %>% 
  mutate(bugs = trimws(bugs, "both")) %>%
  mutate(bugs = gsub("E. cloaceae", "Enterobacter cloaceae", bugs)) %>% 
  mutate(bugs = gsub("Human Adenovirus", "Mastadenovirus", bugs)) %>% 
  separate(bugs, c("bug_genus"), "\\ ", remove = F) %>% 
  mutate(bug_genus = ifelse(grepl("metapneumovirus", bugs), "Metapneumovirus", bug_genus)) %>%
  mutate(bug_genus = ifelse(grepl("rhinovirus", bugs), "Enterovirus", bug_genus)) %>%
  mutate(validity = ifelse(is.na(validity), "Valid", validity)) %>%
  mutate(validity = ifelse(grepl("Partial", validity), "Invalid", validity))

parsed %>%
  fwrite("data/metadata/parsed_microbiology_results.csv")

parsed %>%
  distinct(bugs)
parsed_filt <- parsed %>%
  filter(bugs != "Invalid test") %>%
  filter(bug_genus != "Coliform") %>%
  filter(!grepl("virus", bug_genus)) %>%
  filter(!grepl("virus", bugs)) %>%
  filter(bug_genus != "Pneumocystis") %>%
  # Create field indicating whether species resolution is available
  mutate(species_resolved = !grepl("group|complex|spp|negative|calcoaceticus\\-baumannii", 
                                   bugs, 
                                   ignore.case = T))

parsed_filt %>%
  fwrite("data/metadata/parsed_microbiology_results.bacterial_sp_only.csv")

read_df <- fread("results/tax_classification_out/abundance_matrices/RA.S.zeroed.csv")
taxa_list <- colnames(read_df)
to_keep <- taxa_list[grepl("Enterobacter|calcoaceticus|Proteus|pittii|nosocomialis", taxa_list)]

parsed_filt %>%
  filter(!(bugs %in% c("Negative", "Invalid test"))) %>%
  distinct(bugs) %>% 
  bind_rows(tibble(bugs = to_keep)) %>% 
  filter(!grepl("calcoaceticus-baumannii|Proteus spp.|Enterobacter spp.", bugs)) %>% 
  mutate(bugs = ifelse(grepl("complex|group", bugs), "Enterobacter cloacae complex sp.", bugs)) %>% 
  fwrite("data/metadata/parsed_microbiology_results.species_list.csv")
