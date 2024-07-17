rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)

bug_prefix <- "s_aureus"
bug_prefix <- "e_coli"
# bug_prefix <- "m_catarrhalis"
aln_path <- str_glue("data/trees/gubbins/{bug_prefix}.filtered_polymorphic_sites.fasta")
out_path <- str_glue("data/alignments/{bug_prefix}.core_genes.gubbins_pruned.trimmed.aln")

matrix_to_dnaset <- function(mat) {
  temp_mat <- apply(mat, 1, paste0, collapse = "")
  dnaset <- DNAStringSet(temp_mat, use.names = T)
  return(dnaset)
}

aln <- readDNAStringSet(aln_path)
mat <- as.matrix(aln)

# Mask gappy sites
prop_site_gaps <- apply(mat, 2,
                        function(x) {sum(x %in% c("-", "N")) / nrow(mat)})

sites_to_mask <- prop_site_gaps > 0.20

mat <- mat[, !sites_to_mask]

masked_aln <- matrix_to_dnaset(mat)

writeXStringSet(masked_aln, 
                out_path)

dim(mat)
