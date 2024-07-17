rm(list = ls())
setwd("C:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)
require(ape)

bug_prefix <- "s_aureus"
st_df <- fread(str_glue("results/mlst_out/{bug_prefix}.filt.mlst.tsv")) %>%
  separate(V1, c(rep(NA, 9), "accession"), "/") %>%
  mutate(accession = gsub(".fna", "", accession)) %>%
  dplyr::rename(st_type = V2) %>%
  filter(st_type == 30| accession %in% c("INHALE_FRESH_9-barcode01-bins.1", "INHALE_FRESH_19-barcode06-bins.2",
                                         "INHALE_FRESH_9-barcode02-bins.1"))
st_df
fna <- readDNAStringSet("data/alignments/m_catarrhalis.core_genes.gubbins_pruned.trimmed.aln")
tree <- read.tree("data/trees/gubbins/m_catarrhalis.final_tree.tre")

fna <- readDNAStringSet("data/trees/clonal_frame/s_aureus.core_gene_alignment_filtered.snps.trimmed.vfasttree.recomb_pruned.filtered.fasta")
tree <- read.tree("data/trees/clonal_frame/s_aureus.core_gene_alignment_filtered.snps.trimmed.vfasttree.recomb_pruned.labelled_tree.newick")

meta <- fread("data/metadata/hunt_et_al_v0.2/ena_metadata.tsv")
patient_meta <- fread("data/metadata/parsed_patient_metadata.filt.csv")

parsed_patient <- tibble(sample_accession = names(fna)) %>%
  filter(grepl("Library|INHALE", sample_accession)) %>%
  separate(sample_accession, c("run", "barcode", NA), "\\-", remove = F) %>%
  mutate(run = gsub("\\.no_human|Library|_raw|INHALE_FRESH_|barcode0|barcode", "", run),
         barcode = gsub("\\.no_human|Library|_raw|INHALE_FRESH_|barcode0|barcode", "", barcode)) %>%
  mutate(run = gsub("a|A", "", run, ignore.case = T)) %>%
  mutate(run_id = str_glue("{run}_{barcode}")) %>%
  left_join(patient_meta %>% select(run_id, collection_date = sample_date)) %>%
  mutate(collection_date = as.Date(collection_date, "%d/%m/%Y"))
  
meta_filt <- meta %>% 
  filter(sample_accession %in% names(fna)) %>%
  select(sample_accession, collection_date) %>%
  mutate(collection_date = as.Date(collection_date, "%Y-%m-%d"))

merged <- bind_rows(parsed_patient, meta_filt) %>%
  filter(!is.na(collection_date)) %>%
  filter(sample_accession %in% st_df$accession)

fna_filt <- fna[merged$sample_accession]

to_remove <- deframe(tibble(sample_accession = tree$tip.label) %>%
  filter(!(sample_accession %in% names(fna_filt))))

tree_filt <- drop.tip(tree, to_remove)

meta.match <- tibble(sample_accession = tree_filt$tip.label) %>%
  left_join(merged) %>%
  mutate(collection_date = decimal_date(collection_date)) %>%
  mutate(parsed_labels = str_glue("{sample_accession}|{collection_date}"))

fna_filt <- fna[meta.match$sample_accession]
names(fna_filt) <- meta.match$parsed_labels

tree_filt$tip.label <- meta.match$parsed_labels

# Write
write.tree(tree_filt, "data/alignments/s_aureus.core_genes.clonalframe_pruned.trimmed.date_complete.tree")
writeXStringSet(fna_filt, "data/alignments/s_aureus.core_genes.clonalframe_pruned.trimmed.date_complete.fasta")
merged %>%
  select(sample_accession, collection_date) %>%
  fwrite("data/alignments/s_aureus.core_genes.clonalframe_pruned.trimmed.date_complete.txt")
