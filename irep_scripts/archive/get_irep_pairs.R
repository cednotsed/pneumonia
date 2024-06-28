setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)

high_microbes <- fread("results/qc_out/high_read_samples.csv")
meta <- fread("data/metadata/parsed_patient_metadata.csv")

rank <- "S"

dat <- fread(str_glue("results/metagenomic_out/abundance_matrix.{rank}.tsv")) %>%
  select(-any_of(c("Homo sapiens", "Homo")), -unclassified) %>%
  filter(!grepl("EXTRACTION_TEST", run)) %>%
  as_tibble() %>%
  mutate(run_id = gsub("INHALE_FRESH_|barcode0|barcode", "", run)) %>%
  mutate(run_id = gsub("a|A", "", run_id, ignore.case = T)) %>%
  mutate(run_id = gsub("-", "_", run_id)) %>%
  # Remove runs not in meta
  filter(run_id %in% meta$run_id) %>%
  filter(run_id %in% high_microbes$run_id) %>%
  mutate(run = str_glue("{run}.fastq.gz"))

# Get taxa with min reads
ref_df <- dat %>%
  pivot_longer(!c(run, run_id), names_to = "taxa", values_to = "read_count") %>%
  filter(read_count >= read_threshold) %>%
  filter(!grepl("phage|virus|Nakaseomyces|Candida|Aspergillus|Cryptococcus|Saccharomyces", taxa))

# Get species to taxid mapping
k2_tax <- fread("databases/k2_pluspf_20231009/inspect.txt",
                skip = 7) %>%
  select(taxid = V5, taxa = V6)

taxids <- k2_tax %>%
  filter(taxa %in% ref_df$taxa)

ref_df %>%
  left_join(taxids) %>%
  fwrite("data/metadata/irep/irep_pairs.txt",
         col.names = F,
         eol = "\n")

nrow(taxids) == length(ref_df$taxa)

taxids %>%
  fwrite("data/metadata/irep/references_to_download.txt",
         col.names = F,
         eol = "\n")


taxids
# Get taxid to acc mapping
# k2_map <- fread("databases/k2_pluspf_20231009/seqid2taxid.map.parsed", sep = " ") %>%
#   mutate(V1 = ifelse(grepl("\\|", V1),
#                      str_split(V1, "\\|", simplify = T)[, 3],
#                      V1)) %>%
#   dplyr::rename(acc = V1, taxid = V2)
# 
# 
# test <- k2_map %>%
#   left_join(k2_tax) %>%
#   filter(species %in% ref_df$taxa)
#   # distinct(taxid)
# 
# ref_df$taxa[!(ref_df$taxa %in% test$species)]
# k2_tax %>%
#   filter(species == "Neisseria sp. oral taxon 014")
# 
# k2_map %>%
#   filter(taxid == 641148)
# k2_tax %>%
#   filter(species %in% ref_df$taxa) %>%
#   distinct(taxid)
