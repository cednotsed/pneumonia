rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(doParallel)

blastout <- fread("results/assembly_out/viral_blast_out/all_contigs.renamed.blast.tsv")

colnames(blastout) <- c("qseqid","sseqid","pident","aln_length","mismatch",
                        "gapopen", "qstart", "qend", "sstart", "send",
                        "evalue", "bitscore")

meta <- fread("data/metadata/all_viruses.220723.filt.QCed.csv") %>%
  select(sseqid = accession, family, genus, species)

blastout %>%
  left_join(meta) %>%
  separate(qseqid, c("run_name", "contig"), "\\.") %>% View()
  group_by(run_name) %>%
  distinct(genus) %>%
  View()
