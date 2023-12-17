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
  mutate(run = gsub("INHALE_FRESH_|barcode0|barcode", "", run)) %>%
  mutate(run = gsub("a|A", "", run, ignore.case = T)) %>%
  mutate(run_id = gsub("-", "_", run)) %>%
  select(-run) %>%
  # Remove runs not in meta
  filter(run_id %in% meta$run_id) %>%
  filter(run_id %in% high_microbes$run_id) %>%
  column_to_rownames("run_id") 

# Get average read length
file_dir <- "data/sequencing_summaries/"
file_list <- list.files(file_dir, full.names = T)

read_df <- foreach(file = file_list, .combine = "bind_rows") %do% {
  # file = file_list[1]
  id <- gsub(file_dir, "", file)
  id <- gsub("_sequencing_summary.parsed.csv", "", id)
  
  temp <- fread(file) %>%
    select(sequence_length_template) %>%
    mutate(run = id)
  
  return(temp)
}

read_df %>%
  summarise(n_runs = n_distinct(run))

median_read_length <- deframe(read_df %>%
  summarise(median_length = median(sequence_length_template)))

median_read_length

# Get average genome_size
patric <- fread("../ARG_host_jumps/data/metadata/genome_metadata/BVBRC_genome.221023.csv") %>%
  rename_all(~tolower(gsub(" ", "_", .x)))

median_genome_size <- deframe(patric %>%
  filter(genome_quality == "Good") %>%
  summarise(median_genome_size = median(size)))

# Calculate min reads
cov_threshold <- 1
read_threshold <- signif(1 * median_genome_size / median_read_length, 2)

# Get taxa with min reads
ref_df <- dat %>%
  pivot_longer(everything(), names_to = "taxa", values_to = "read_count") %>%
  filter(read_count >= read_threshold) %>%
  distinct(taxa) %>%
  filter(!grepl("phage|virus|Nakaseomyces|Candida|Aspergillus|Cryptococcus|Saccharomyces", taxa))

# Get species to taxid mapping
k2_tax <- fread("databases/k2_pluspf_20231009/inspect.txt",
                skip = 7) %>%
  select(taxid = V5, species = V6)

length(ref_df$taxa)

k2_tax %>%
  filter(grepl("P"))
taxids <- k2_tax %>%
  filter(species %in% ref_df$taxa) %>% View()
  distinct(taxid)

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
