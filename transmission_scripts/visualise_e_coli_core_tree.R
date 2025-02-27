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
require(phytools)

# Get tree metadata
meta <- fread("data/metadata/parsed_patient_metadata.filt.csv")

bug_prefix <- "e_coli"
tree <- read.tree(str_glue("data/trees/gubbins/{bug_prefix}.final_tree.tre"))
rooted <- midpoint.root(tree)

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
  head(15)

st_to_keep <- st_df %>%
  filter(grepl("Library|INHALE", accession)) %>%
  bind_rows(st_filt) %>%
  distinct(st_type)

# Match metadata to tips
meta.match <- tibble(accession = rooted$tip.label) %>%
  separate(accession, c("run", "barcode"), "-", remove = F) %>%
  mutate(run = gsub("Library|INHALE_FRESH_", "", run),
         barcode = gsub("barcode|barcode0", "", barcode)) %>%
  mutate(barcode = gsub("a", "", barcode)) %>%
  mutate(run_id = str_glue("{run}_{barcode}")) %>%
  left_join(meta) %>%
  left_join(st_df) %>%
  mutate(is_assembly = !is.na(hap_vap_cap)) %>%
  mutate(annot = str_glue("{accession}|{hospital}|{sample_date}|{hap_vap2}|ST_TYPE={st_type}")) %>%
  mutate(sample_date = as.Date(sample_date, "%d/%m/%Y")) %>%
  mutate(st_annot = ifelse(st_type %in% st_to_keep$st_type, st_type, NA)) %>%
  mutate(st_annot = ifelse(st_type == "-", NA, st_annot)) %>%
  mutate(st_text = ifelse(is_assembly, st_annot, NA))

annot_tree <- rooted
annot_tree$tip.label <- meta.match$annot
write.tree(annot_tree, str_glue("data/trees/gubbins/{bug_prefix}.gubbins.annot.newick"))

all(tree$tip.label == meta.match$accession)

# Remove non-HAP VAP isolates
new <- rooted$tip.label[grepl("INHALE|Lib", tree$tip.label)]
to_remove <- c("SAMEA8156163", "SAMEA8156220")
# # to_remove <- meta.match %>%
#   # filter(!()))
# 
rooted_filt <- drop.tip(rooted, to_remove)
meta_filt <- meta.match %>%
  filter(!(accession %in% to_remove))

pal <- list(UCLH = "darkcyan", RFH = "tan1", Cromwell = "darkorchid4", 
            GSTT = "#2e4057", `Local GP (controls)` = "khaki", GOSH = "cornflowerblue",
            `North Mid` = "#d1495b", `ChelWest` = "cadetblue2")


st_pal <- distinctColorPalette(n_distinct(meta.match$st_annot))

width <- 7000

ggtree(rooted_filt,
       size = 0.1,
       color = "grey40",
       ladderize = T,
       options(ignore.negative.edge = TRUE)) %<+% meta_filt +
  geom_tippoint(aes(fill = hospital, color = ifelse(!is.na(hospital), "black", NA)), 
                size = 3,
                alpha = 0.9,
                pch = 21) +
  geom_tiplab(aes(label = st_text), hjust = 0) +
  scale_color_manual(values = c("black", NA), na.value = NA, guide = 'none') +
  scale_fill_manual(values = pal, na.translate = F) +
  labs(fill = "Site") +
  new_scale_fill() +
  geom_fruit(geom = geom_tile,
             aes(fill = st_annot),
             offset = 0.08,
             width = width) +
  scale_fill_manual(values = st_pal, na.translate = F) +
  labs(fill = "MLST") +
  theme_tree2() +
  coord_cartesian(clip="off")

ggsave(str_glue("results/transmission_out/{bug_prefix}.recomb_pruned_tree.pdf"), width = 4, height = 4, dpi = 300)


