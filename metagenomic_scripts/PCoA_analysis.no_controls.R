setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(vegan)
require(ape)

meta <- fread("data/metadata/parsed_patient_metadata.csv")

rank <- "S"
RA_filt <- fread(str_glue("results/metagenomic_out/RA.{rank}.zeroed.csv")) %>%
  left_join(meta %>% select(run_id, hap_vap_cap)) %>%
  filter(hap_vap_cap != "Water control") %>%
  select(-hap_vap_cap) %>%
  column_to_rownames("run_id")

# Bray-Curtis PCoA
bc <- vegdist(as.matrix(RA_filt), method = "bray")

PCOA <- pcoa(bc)

eigenvals <- round(PCOA$values$Relative_eig * 100, 1)

plot_df <- as.data.frame(PCOA$vectors) %>%
  rownames_to_column("run_id") %>%
  left_join(meta)

plot_df %>%
  ggplot(aes(Axis.1, Axis.2, color = hap_vap_cap, shape = hap_vap_cap)) +
  geom_point() +
  scale_color_manual(values = c("indianred", "blue", "goldenrod", "grey50")) +
  theme_classic() + 
  labs(x = str_glue("PCo1 ({eigenvals[1]}%)"),
       y = str_glue("PCo2 ({eigenvals[2]}%)"),
       color = "Sample Type",
       shape = "Sample Type",
       title = "Bray-Curtis PCoA")

# ggsave(str_glue("results/metagenomic_out/PCoA.{rank}.no_decontamination.png"),
#        height = 5, width = 5)

# pca <- prcomp(RA_filt %>%
#                 rownames_to_column("run_id") %>%
#                 filter(!(run_id %in% c("31_4", "22_8"))) %>%
#                 column_to_rownames("run_id"), 
#               retx = T)
# 
# pca_eigenvals <- round(pca$sdev^2 / sum(pca$sdev^2) * 100, 1)
# 
# as.data.frame(pca$x) %>%
#   rownames_to_column("run_id") %>%
#   left_join(meta) %>%
#   ggplot(aes(x = PC1, y = PC3, color = hap_vap_cap, shape = hap_vap_cap)) +
#   geom_point() +
#   # geom_text(aes(label = run_id)) +
#   scale_color_manual(values = c("indianred", "blue", "goldenrod", "grey50"))+
#   theme_classic() + 
#   labs(x = str_glue("PCo1 ({pca_eigenvals[1]}%)"),
#        y = str_glue("PCo2 ({pca_eigenvals[2]}%)"),
#        color = "Sample Type",
#        shape = "Sample Type",
#        title = "PCA")
#  
# long_df <- RA_filt %>%
#   rownames_to_column("run_id") %>%
#   pivot_longer(!run_id, names_to = "genus", values_to = "rel_a") %>%
#   left_join(meta)
# 
# summary_stats <- long_df %>%
#   group_by(genus) %>%
#   summarise(max_abundance = max(rel_a)) %>%
#   arrange(desc(max_abundance)) %>%
#   head(30)
# 
# pal <- distinctColorPalette(n_distinct(summary_stats$genus))
# long_df %>%
#   filter(genus %in% summary_stats$genus) %>%
#   mutate(genus = factor(genus, unique(summary_stats$genus))) %>%
#   ggplot(aes(x = run_id, y = rel_a, fill = genus)) +
#   geom_bar(stat = "identity", 
#            position = "stack",
#            color = "black") +
#   facet_grid(~hap_vap_cap, scales = "free", space = "free") +
#   theme_classic() +
#   scale_fill_manual(values = pal) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# # dat
# # RA_parsed[1:5, 1:5]
# # RA_long <- RA_parsed %>%
# #   pivot_longer(!c("run", "barcode"), 
# #                names_to = "species",
# #                values_to = "rel_a")
# # merged_df <- meta %>%
# #   left_join(RA_long)
# #   
# # merged_df %>% 
# #   filter(micro_organism == "Morganella morganii") %>%
# #   filter(grepl("Morganella", species))
# # 
# # RA_df %>%
# #   rownames_to_column("run") %>%
# #   filter(run == "INHALE_FRESH_1-barcode01") %>%
# #   select(`Morganella morganii`)
# #   filter()
# #   
# #   
# # meta %>%
# #   distinct(curetis_valid)
# # dat %>%
# #   pivot_longer(!c("run", "barcode"), 
# #                names_to = "species", 
# #                values_to = "read_count") %>%
# #   filter(species == "Homo sapiens") %>%
# #   arrange(desc(read_count))
# # meta %>%
# #   select(run, barcode, )
# # meta %>%
# #   distinct(run)
# # colnames(meta)
# 
# # # Parse IDs
# # RA_parsed <- RA_filt %>%
# #   # Remove IDs with zero reads
# #   filter(is.na(`Pseudomonas sp. T1-3-2`))
# #   # rownames_to_column("run_id")
# #   # separate(run, into = c("run", "barcode"), "-") %>%
# #   # mutate(barcode = as.numeric(gsub("barcode", "", barcode))) %>%
# #   # mutate(barcode = as.character(barcode)) %>%
# #   # mutate(barcode = replace_na(barcode, "12a")) %>%
# #   # mutate(run = gsub("INHALE_FRESH_", "", run)) %>%
# #   # mutate(run_id = str_glue("{run}_{barcode}")) %>%
# #   # select(-run, -barcode) %>%
# #   # relocate(run_id, .before = 1) %>%
# #   # column_to_rownames("run_id")