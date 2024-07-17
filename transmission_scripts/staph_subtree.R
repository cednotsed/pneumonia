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

bug_prefix <- "s_aureus"
tree <- read.tree(str_glue("data/trees/clonal_frame/{bug_prefix}.core_gene_alignment_filtered.snps.trimmed.vfasttree.recomb_pruned.labelled_tree.newick"))

# ST types 
st_df <- fread(str_glue("results/mlst_out/{bug_prefix}.filt.mlst.tsv")) %>%
  separate(V1, c(rep(NA, 9), "accession"), "/") %>%
  mutate(accession = gsub(".fna", "", accession)) %>%
  dplyr::rename(st_type = V2)

# Get top 10 ST types
st_filt <- st_df %>%
  group_by(st_type) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>% 
  head(10)

st_to_keep <- st_df %>%
  filter(grepl("Library|INHALE", accession)) %>%
  bind_rows(st_filt) %>%
  distinct(st_type)

# Match metadata to tips
meta.match <- tibble(accession = tree$tip.label) %>%
  separate(accession, c("run", "barcode"), "-", remove = F) %>%
  mutate(run = gsub("Library|INHALE_FRESH_", "", run),
         barcode = gsub("barcode|barcode0", "", barcode)) %>%
  mutate(barcode = gsub("a", "", barcode)) %>%
  mutate(run_id = str_glue("{run}_{barcode}")) %>%
  left_join(meta) %>%
  left_join(st_df) %>%
  mutate(is_assembly = !is.na(hap_vap_cap)) %>%
  mutate(annot = str_glue("{accession}|{hospital}|{sample_date}|{hap_vap_cap}|ST_TYPE={st_type}")) %>%
  mutate(sample_date = as.Date(sample_date, "%d/%m/%Y")) %>%
  mutate(st_annot = ifelse(st_type %in% st_to_keep$st_type, st_type, NA)) %>%
  mutate(st_annot = ifelse(st_type == "-", NA, st_annot)) %>%
  mutate(st_text = ifelse(is_assembly, st_annot, NA))

of_interest <- fread("data/trees/clonal_frame/s_aureus_interest.txt", header = F)$V1

to_remove <- meta.match$accession[!(meta.match$accession %in% of_interest)]

tree_filt <- drop.tip(tree, to_remove)

meta_filt <- tibble(accession = tree_filt$tip.label) %>%
  left_join(meta.match)

pal <- list(UCLH = "darkcyan", RFH = "tan1", Cromwell = "darkorchid4", 
            GSTT = "#2e4057", `Local GP (controls)` = "khaki", GOSH = "cornflowerblue",
            `North Mid` = "#d1495b", `ChelWest` = "cadetblue2")


st_pal <- distinctColorPalette(n_distinct(meta.match$st_annot))

ggtree(tree_filt,
       size = 0.1,
       # layout = "fan",
       branch.length = "none",
       color = "grey40",
       options(ignore.negative.edge = TRUE)) %<+% meta.match +
  geom_tippoint(aes(fill = hospital, 
                    shape = hap_vap_cap), 
                size = 3,
                alpha = 0.9) +
  geom_tiplab(aes(label = sample_date)) +
  scale_shape_manual(values = c(21, 24)) +
  scale_color_manual(values = c("black", NA), 
                     na.value = NA, guide = 'none') +
  scale_fill_manual(values = pal) +
  labs(fill = "Site")

ggsave(str_glue("results/transmission_out/{bug_prefix}.recomb_pruned_tree.pdf"), width = 4, height = 4, dpi = 300)

# # Plot tree
# col_pal <- c("darkcyan", "seagreen4", "darkgoldenrod3", 
#              "#2e4057", "#d1495b", "cornflowerblue", 
#              "maroon4", "darkorchid4", "slateblue4",
#              "black", "darkgrey", "lightpink", 
#              "navajowhite3", "darkolivegreen4", "lightgreen",
#              "tan2", "blue", "firebrick4")
# 
# set.seed(66)
# random_pal <- distinctColorPalette(length(unique(meta.match$genus)))
# # random_pal[length(random_pal) - 1] <- "darkolivegreen4"
# 
# ggtree(rooted,
#        size = 0.001,
#        # branch.length = "none",
#        color = "darkslategrey",
#        options(ignore.negative.edge = TRUE)) %<+% meta.match +
#   geom_tippoint(aes(hjust = 0.5, color = genus), 
#                 alpha = 1, 
#                 size = 1) +
#   scale_color_discrete(na.translate = F) +
#   geom_fruit(geom = geom_tile,
#              aes(fill = host_order),
#              offset = 0.15,
#              width = 0.03) +
#   labs(color = "Viral genus", 
#        fill = "Host order") +
#   scale_fill_manual(values = col_pal,
#                     na.translate = F)
# 
# ggsave("results/phylogenetic_out/cov_NJ_mash.rooted.linear.pdf",
#        dpi = 300,
#        width = 12,
#        height = 12)
# 
# ggtree(rooted,
#        size = 0.001,
#        # branch.length = "none",
#        color = "darkslategrey",
#        options(ignore.negative.edge = TRUE)) %<+% meta.match +
#   geom_tippoint(aes(hjust = 0.5, color = genus), 
#                 alpha = 1, 
#                 size = 1) +
#   scale_color_discrete(na.translate = F) +
#   geom_tiplab(aes(label = tip_label),
#               size = 0.5)
# ggsave("results/phylogenetic_out/cov_NJ_mash.rooted.linear.tip_labels.pdf",
#        dpi = 100,
#        width = 30,
#        height = 30)
