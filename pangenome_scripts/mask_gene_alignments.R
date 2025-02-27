rm(list = ls())
setwd("/mnt/c/git_repos/pneumonia/")
require(tidyverse, quietly = T)
# require(data.table, quietly = T)
# require(foreach, quietly = T)
require(Biostrings, quietly = T)

matrix_to_dnaset <- function(mat) {
  temp_mat <- apply(mat, 1, paste0, collapse = "")
  dnaset <- DNAStringSet(temp_mat, use.names = T)
  return(dnaset)
}

mask_aln <- function(aln) {
  mat <- as.matrix(aln)
  
  # Mask gappy sites
  prop_site_gaps <- apply(mat, 2,
                          function(x) {sum(x %in% c("-", "N")) / nrow(mat)})
  
  site_to_mask <- prop_site_gaps > 0.20
  
  mat <- mat[, !site_to_mask]
  
  masked_aln <- matrix_to_dnaset(mat)
  
  print(sum(site_to_mask))
  return(masked_aln)  
}



# Run
args = commandArgs(trailingOnly=TRUE)
bug_prefix <- args[1]
idx <- as.numeric(args[2])

# Read full snp alignment
snp_aln <- readDNAStringSet(str_glue("data/alignments/{bug_prefix}.core_gene_alignment_filtered.snps.trimmed.aln"))
to_keep <- names(snp_aln)

# Read gene list
file_dir <- str_glue("results/pangenome_out/{bug_prefix}.filt/aligned_gene_sequences/")
file_list <- list.files(file_dir, full.names = T)
file_name <- file_list[idx]

out_dir <- str_glue("results/pangenome_out/{bug_prefix}.filt/aligned_gene_sequences.masked")
dir.create(out_dir)

# MAIN
prefix <- gsub(file_dir, "", file_name)
out_path <- gsub(".fas", ".masked.fna", prefix)

gene_aln <- readDNAStringSet(file_name)
names(gene_aln) <- str_split(names(gene_aln), "\\;", simplify = T)[, 1]

# Create blanks for missing genes
missing <- names(snp_aln)[!(names(snp_aln) %in% names(gene_aln))]
missing_mat <- matrix("-", nrow = length(missing), ncol = unique(width(gene_aln)))
missing_aln <- matrix_to_dnaset(missing_mat)
names(missing_aln) <- missing

# Filter gene alignment
aln_filt <- c(gene_aln, missing_aln)[names(snp_aln)]
all(names(aln_filt) == names(snp_aln))

# Mask gene alignment
masked_aln <- mask_aln(aln_filt)

writeXStringSet(masked_aln, str_glue("{out_dir}/{out_path}"))
