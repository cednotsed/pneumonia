rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(ape)
require(randomcoloR)
require(lubridate)

# Get tree metadata
sp_df <- fread("results/assembly_out/assembly_merged_metadata.csv")
meta <- fread("data/metadata/parsed_patient_metadata.filt.csv") %>%
  filter(hap_vap_cap %in% c("HAP", "VAP"))

bug_prefix <- "s_aureus"
bug_prefix <- "e_coli"

staph_tree <- read.tree(str_glue("data/trees/clonal_frame/s_aureus.core_gene_alignment_filtered.snps.trimmed.vfasttree.recomb_pruned.labelled_tree.newick"))
ecoli_tree <- read.tree(str_glue("data/trees/clonal_frame/e_coli.core_gene_alignment_filtered.snps.trimmed.vfasttree.recomb_pruned.labelled_tree.newick"))
hinfluenzae_tree <- read.tree(str_glue("data/trees/clonal_frame/h_influenzae.core_gene_alignment_filtered.snps.trimmed.vfasttree.recomb_pruned.labelled_tree.newick"))

# ST types 
st_df <- bind_rows(fread(str_glue("results/mlst_out/s_aureus.filt.mlst.tsv")),
                   fread(str_glue("results/mlst_out/e_coli.filt.mlst.tsv")),
                   fread(str_glue("results/mlst_out/h_influenzae.filt.mlst.tsv"))) %>%
  separate(V1, c(rep(NA, 9), "accession"), "/") %>%
  mutate(accession = gsub(".fna", "", accession)) %>%
  dplyr::rename(st_type = V2)

# Match metadata to tips
meta.match <- tibble(accession = c(staph_tree$tip.label, ecoli_tree$tip.label, hinfluenzae_tree$tip.label)) %>%
  separate(accession, c("run", "barcode"), "-", remove = F) %>%
  mutate(run = gsub("Library|INHALE_FRESH_", "", run),
         barcode = gsub("barcode|barcode0", "", barcode)) %>%
  mutate(barcode = gsub("a", "", barcode)) %>%
  mutate(run_id = str_glue("{run}_{barcode}")) %>%
  filter(!grepl("SAM", run_id)) %>%
  left_join(meta) %>% 
  left_join(st_df) %>%
  left_join(sp_df %>% select(accession = genome_name, species)) %>%
  mutate(species = gsub("_D|_E", "", species)) %>%
  mutate(sample_date = as.Date(sample_date, "%d/%m/%Y")) %>%
  filter(!is.na(hospital))

pal <- list(UCLH = "darkcyan", RFH = "tan1", Cromwell = "darkorchid4", 
            GSTT = "#2e4057", `Local GP (controls)` = "khaki", GOSH = "cornflowerblue",
            `North Mid` = "#d1495b", `ChelWest` = "cadetblue2")

meta.match %>%
  ggplot(aes(x = sample_date, y = hospital, fill = hospital, shape = hap_vap2)) +
  geom_point(color = "black",
             size = 3,
             position = position_jitter(width = 0, height = 0.3)) +
  scale_shape_manual(values = c(21, 22, 24)) +
  geom_text(aes(label = st_type)) +
  scale_fill_manual(values = pal, guide = "none") + 
  facet_grid(rows = vars(species)) +
  theme_bw() +
  scale_x_date(breaks = "1 month", date_labels = "%b-%y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.y = element_blank()) +
  labs(x = "Sample collection date", y = "Site")

ggsave("results/transmission_out/date_distribution.pdf", dpi = 600, width = 8, height = 4)

# comb_df <- t(combn(novel$accession, 2))
# mat <- cophenetic.phylo(tree)
# 
# morsels <- foreach(i = seq(nrow(comb_df))) %do% {
#   # i = 1
#   acc1 <- comb_df[i, 1]
#   acc2 <- comb_df[i, 2]
#   
#   cophen_dist <- mat[acc1, acc2]
#   date_dist <- abs(difftime((novel %>% 
#                           filter(accession == acc1))$sample_date,
#                         (novel %>% 
#                            filter(accession == acc2))$sample_date))
#   
#   tibble(cophen_dist = cophen_dist, 
#          time_diff = date_dist)
#   
# }
# 
# plot_df <- bind_rows(morsels)
# 
# cor.test(plot_df$cophen_dist, as.numeric(plot_df$time_diff), method = "spearman")
# 
# plot_df %>%
#   ggplot(aes(x = time_diff, y = cophen_dist)) +
#   geom_point()
