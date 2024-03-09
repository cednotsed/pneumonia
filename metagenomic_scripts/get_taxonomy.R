setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(vegan)
require(hillR)
require(ggpubr)

meta <- fread("data/metadata/parsed_patient_metadata.csv")

dat <- fread(str_glue("results/metagenomic_out/abundance_matrix.{rank}.tsv")) %>%
  # select(-`Homo sapiens`, -unclassified) %>%
  filter(!grepl("EXTRACTION_TEST", run)) %>%
  as_tibble() %>%
  mutate(run = gsub("INHALE_FRESH_|barcode0|barcode", "", run)) %>%
  mutate(run_id = gsub("-", "_", run)) %>%
  select(-run) %>%
  # Remove runs not in meta
  filter(run_id %in% meta$run_id) %>%
  column_to_rownames("run_id") 

taxa_list <- fread("databases/k2_pluspf_20231009/inspect.txt",
                   skip = 7)

taxa_filt <- taxa_list %>% 
  filter(V4 %in% c("D", "K", "P", "S")) 

morsels <- foreach(i = seq(nrow(taxa_filt))) %do% {
  # i = 1
  row <- taxa_filt[i, ]
  if(row$V4 == "D") {
    domain_name <- row$V6
    return(NULL)
  } else {
    if(row$V4 == "P") {
      phylum_name <- row$V6
      return(NULL)
    } else {
      return(row %>% mutate(domain = domain_name,
                            phylum = phylum_name))
    }
  }
}

domain_df <- bind_rows(morsels) %>%
  mutate(group = case_when(domain == "Bacteria" ~ "Bacteria",
                           domain == "Archaea" ~ "Archaea",
                           domain == "Viruses" ~ "Viruses",
                           domain == "Eukaryota" &
                             phylum == "Chordata" ~ "Humans",
                           domain == "Eukaryota" &
                             phylum %in% c("Ascomycota", "Basidiomycota") ~ "Fungi",
                           TRUE ~ "Other eukaryotes"))

domain_filt <- domain_df %>%
  select(species = V6, domain, phylum, group)

fwrite(domain_filt, "data/metadata/k2_database_taxonomy.csv")

domain_filt
# domain_df %>% 
#   filter(group == "Fungi") %>% View()
# eukaryote_df <- domain_df %>%
#   filter(domain == "Eukaryota")
# 
# 
# eukaryote_morsels <- foreach(i = seq(nrow(eukaryote_df))) %do% {
#   # i = 1
#   row <- eukaryote_df[i, ]
#   if(row$V4 == "K") {
#     kingdom_name <- row$V6
#     return(NULL)
#   } else {
#     return(row %>% mutate(kingdom = kingdom_name))
#   }
# }
# 
# bind_rows(eukaryote_morsels) %>%
#   distinct(domain, kingdom, phylum) %>%
#   View()
# 
# domain_df %>%
#   filter(domain == "Eukaryota") %>%
#   distinct(phylum)
#   filter(grepl("Plasmodium", V6))
#   filter(!(domain %in% c("Archaea", "Bacteria", "Viruses"))) %>%
#   distinct(domain)
#   distinct(phylum)
#   filter(domain == "Eukaryota") %>% View()
# taxa_list %>% 
#   filter(V4 == "K")
# 
# taxa_list %>% 
#   filter(grepl("Plasmodium", V6))
