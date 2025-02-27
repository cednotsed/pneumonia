rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(ape)
require(ggtree)
require(ggtreeExtra)
require(ggnewscale)
require(ggutils)
require(randomcoloR)

# Get tree metadata
meta <- fread("data/metadata/parsed_patient_metadata.filt.csv")

to_remove <- deframe(meta %>%
  filter(hap_vap_cap %in% c("CAP", "Healthy")) %>%
  distinct(run_id))

bug_prefix <- "m_catarrhalis"
bug_prefix <- "s_aureus"
bug_prefix <- "e_coli"

morsels <- foreach(bug_prefix = c("m_catarrhalis", "s_aureus", "e_coli")) %do% {
  tree <- read.tree(str_glue("data/trees/clonal_frame/{bug_prefix}.core_gene_alignment_filtered.snps.trimmed.vfasttree.recomb_pruned.labelled_tree.newick"))
  
  # Match metadata to tips
  meta.match <- tibble(accession = tree$tip.label) %>%
    separate(accession, c("run", "barcode"), "-", remove = F) %>%
    mutate(run = gsub("Library|INHALE_FRESH_", "", run),
           barcode = gsub("barcode|barcode0", "", barcode)) %>%
    mutate(barcode = gsub("a", "", barcode)) %>%
    mutate(run_id = str_glue("{run}_{barcode}")) %>%
    left_join(meta) %>%
    mutate(is_assembly = !is.na(hap_vap_cap))
  
  accession_to_remove <- deframe(meta.match %>%
      filter(run_id %in% to_remove) %>%
      select(accession))
  
  # Remove non HAP or VAP
  tree_filt <- drop.tip(tree, accession_to_remove)
  
  # Get cophenetic distance
  cophen_mat <- cophenetic.phylo(tree_filt)
  
  # Novel genomes
  novel_names <- deframe(meta.match %>%
      filter(is_assembly) %>%
      filter(!(accession %in% accession_to_remove)) %>%
      select(accession))
  
  mat_filt <- cophen_mat[novel_names, novel_names]
  obs_dist <- mean(mat_filt[upper.tri(mat_filt, diag = F)])
  
  # Permutation
  all_names <- colnames(cophen_mat)
  
  crumbs <- foreach(i = seq(1000)) %do% {
    temp_names <- sample(all_names, length(novel_names), replace = F) 
    temp_filt <- cophen_mat[temp_names, temp_names]
    temp_dist <- mean(temp_filt[upper.tri(temp_filt, diag = F)])
    
    return(tibble(distance = temp_dist))
  }

  perm_results <- bind_rows(crumbs) %>%
    mutate(bug_prefix = bug_prefix) %>%
    mutate(obs_dist = obs_dist) 
  
  return(perm_results)
}

merged <- bind_rows(morsels)

obs_df <- merged %>% distinct(bug_prefix, obs_dist)
p_df <- merged %>%
  group_by(bug_prefix) %>%
  summarise(n = sum(distance <= obs_dist) / 1000)
merged %>%  
  ggplot(aes(x = bug_prefix, y = distance)) +
  geom_violin() +
  geom_point(aes(x = bug_prefix, y = obs_dist),
             data = obs_df)
