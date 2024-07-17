rm(list = ls())
setwd("C:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)

dom_df <- fread("results/benchmarking_out/undiagnosed_dominant_bugs.csv")
  # dplyr::rename(dom_species = species)
dom_df
tax_meta <- fread("databases/k2_pluspfp_20240112/inspect.txt") %>%
  select(taxon = V6, taxid = V5) %>%
  mutate(taxid = as.character(taxid))

irep_df <- fread("results/irep_out/irep_results.tsv") %>%
  select(run, Ori, Ter, bPTR) %>%
  separate(run, c("run", "taxid"), "\\.") %>%
  mutate(run = gsub("\\.no_human|Library|_raw|INHALE_FRESH_|barcode0|barcode", "", run)) %>%
  mutate(run = gsub("a|A", "", run, ignore.case = T)) %>%
  mutate(run_id = gsub("-", "_", run)) %>%
  select(-run) %>%
  left_join(tax_meta) %>%
  filter(bPTR != "n/a") %>%
  dplyr::rename(species = taxon)

irep_df

dom_filt <- dom_df %>%
  filter(!grepl("Candida|Aspergillus|Penicillium", species))

dom_filt %>%
  left_join(irep_df) %>% View()
  left_join(df %>%
              left_join(tax_meta)) %>% View()
  filter(bPTR != "NA") %>% View()

  dom_df$species[!(dom_df$species %in% irep_df$species)]
  