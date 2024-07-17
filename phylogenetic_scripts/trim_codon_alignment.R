rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)

# bug_prefix <- "s_aureus"
# bug_prefix <- "e_coli"
bug_prefix <- "m_catarrhalis"
aln_path <- str_glue("results/pangenome_out/{bug_prefix}.filt/{bug_prefix}.core_gene_alignment_filtered.snps.aln")
out_path <- str_glue("data/alignments/{bug_prefix}.core_gene_alignment_filtered.snps.trimmed.aln")

matrix_to_dnaset <- function(mat) {
  temp_mat <- apply(mat, 1, paste0, collapse = "")
  dnaset <- DNAStringSet(temp_mat, use.names = T)
  return(dnaset)
}

aln <- readDNAStringSet(aln_path)
mat <- as.matrix(aln)

# Remove gappy sequences
prop_gaps <- apply(mat, 1,
                   function(x) {sum(x %in% c("-", "N")) / ncol(mat)})

to_remove <- prop_gaps > 0.20
prop_gaps[to_remove]
print(names(aln)[to_remove])
mat <- mat[!to_remove, ]

dim(mat)

# Mask gappy sites
prop_site_gaps <- apply(mat, 2,
                   function(x) {sum(x %in% c("-", "N")) / nrow(mat)})

site_to_mask <- prop_site_gaps > 0.20

mat <- mat[, !site_to_mask]

masked_aln <- matrix_to_dnaset(mat)

writeXStringSet(masked_aln, 
                out_path)
