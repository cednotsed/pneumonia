rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)

ref_df <- fread("data/metadata/irep/undiagnosed_bugs.csv")

# Get species to taxid mapping
k2_tax <- fread("databases/k2_pluspf_20231009/inspect.txt",
                skip = 7) %>%
  select(taxid = V5, species = V6)

taxids <- k2_tax %>%
  filter(species %in% ref_df$species)

taxids %>%
  select(taxid) %>%
  fwrite("data/metadata/irep/taxids_to_download.txt",
         col.names = F,
         eol = "\n")

# After downloading references get reference paths
ref_df %>%
  distinct(run_id) %>%
  fwrite("data/metadata/irep/runs_to_transfer.txt",
         eol = "\n",
         col.names = F)

ref_list <- list.files("data/genomes/irep_references/", recursive = T)
ref_meta <- as.data.frame(str_split(ref_list, "\\/", simplify = T)) %>%
  mutate(ref_path = str_glue("/mnt/c/git_repos/pneumonia/data/genomes/irep_references/{ref_list}")) %>%
  select(taxid = V1, ref_path) %>%
  mutate(taxid = as.numeric(taxid))

# Save pairs
ref_df %>%
  left_join(k2_tax) %>%
  select(-read_count) %>%
  left_join(ref_meta) %>%
  separate(run_id, c("run", "index"), "\\_", remove = F) %>%
  mutate(index = as.numeric(index)) %>%
  mutate(index = ifelse(index < 10, str_glue("0{index}"), as.character(index))) %>%
  mutate(index = ifelse(index == 12, str_glue("{index}a"), index)) %>%
  mutate(ref_path = gsub(".fna.gz", ".filt.fna", ref_path)) %>%
  mutate(fq_path = str_glue("/mnt/c/git_repos/pneumonia/data/basecalled_fastqs/no_humans/INHALE_FRESH_{run}-barcode{index}.no_human.fastq.gz")) %>% 
  select(-run, -index) %>% 
  fwrite("data/metadata/irep/irep.config.tsv",
         eol = "\n",
         sep = "\t",
         col.names = F)

file_list <- list.files("data/genomes/irep_references/", recursive = T, full.names = T)
morsels <- foreach(path = file_list) %do% {
  fna <- readDNAStringSet(path)
  
  n_plasmids <- sum(grepl("plasmid", names(fna)))
  n_chrom <- sum(grepl("chromosome", names(fna)))
  n_acc <- length(fna)
  
  fna_filt <- fna[1]
  writeXStringSet(fna_filt, gsub(".fna", ".filt.fna", path))
  # tibble(path = path, n_acc = n_acc, n_chrom = n_chrom, n_plasmids = n_plasmids)
}

bind_rows(morsels) %>% View()
