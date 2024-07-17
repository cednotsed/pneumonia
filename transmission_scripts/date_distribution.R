rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(ape)
require(randomcoloR)
require(lubridate)

# Get tree metadata
meta <- fread("data/metadata/parsed_patient_metadata.filt.csv")

bug_prefix <- "m_catarrhalis"
bug_prefix <- "s_aureus"
bug_prefix <- "e_coli"

tree <- read.tree(str_glue("data/trees/clonal_frame/{bug_prefix}.core_gene_alignment_filtered.snps.trimmed.vfasttree.recomb_pruned.labelled_tree.newick"))

# M catarrhalis to remove
to_remove <- fread("data/trees/m_catarrhalis.to_remove.txt", header = F)$V1
tree <- drop.tip(tree, to_remove)

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

novel <- meta.match %>%
  filter(!is.na(hap_vap_cap)) %>%
  select(accession, hap_vap_cap, sample_date, st_type, hospital)

pal <- list(UCLH = "darkcyan", RFH = "tan1", Cromwell = "darkorchid4", 
            GSTT = "#2e4057", `Local GP (controls)` = "khaki", GOSH = "cornflowerblue",
            `North Mid` = "#d1495b", `ChelWest` = "cadetblue2")

novel %>%
  ggplot(aes(x = sample_date, y = hospital, fill = hospital, shape = hap_vap_cap)) +
  geom_point(color = "black",
             size = 3,
             position = position_jitter(width = 0, height = 0.3)) +
  scale_shape_manual(values = c(21, 22, 24)) +
  geom_text(aes(label = st_type)) +
  scale_fill_manual(values = pal, guide = "none") +
  theme_bw() +
  scale_x_date(breaks = "1 month", date_labels = "%b-%y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.y = element_blank()) +
  labs(x = "Sample collection date", y = "Site")

ggsave("results/transmission_out/date_distribution.pdf", dpi = 600, width = 5, height = 3)

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
