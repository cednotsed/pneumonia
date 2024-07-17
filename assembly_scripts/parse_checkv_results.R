rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(doParallel)
require(ape)

# CheckV results
check_df <- fread("results/assembly_out/checkv_out/quality_summary.tsv")

check_filt <- check_df %>%
  filter(completeness >= 90, 
         contamination <= 5) %>% 
  filter(warnings == "") %>% 
  filter(provirus == "No")

# Parse blast results
blastout <- fread("results/assembly_out/viral_blast_out/all_contigs.renamed.blast.tsv")

colnames(blastout) <- c("qseqid","sseqid","pident","aln_length","mismatch",
                        "gapopen", "qstart", "qend", "sstart", "send",
                        "evalue", "bitscore")

meta <- fread("data/metadata/all_viruses.220723.filt.QCed.csv") %>%
  select(sseqid = accession, family, genus, species)

merged <- blastout %>%
  left_join(meta)

merged %>% 
  filter(family == "Anelloviridae") %>%
  arrange(desc(aln_length))
View(merged)
  separate(qseqid, c("run_name", "contig"), "\\.", remove = F) %>%
  dplyr::rename(contig_id = qseqid) %>%
  right_join(check_filt)

merged$species
group_by(run_name) %>%
  distinct(genus) %>%
  View()

check_filt %>%
  filter(contig_id %in% unique(blastout$qseqid))

unique(blastout$qseqid) 
unique(check_filt$contig_id)
