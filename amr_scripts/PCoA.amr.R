rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(vegan)
require(ape)

meta_filt <- fread("data/metadata/parsed_patient_metadata.filt.csv") %>%
  filter(high_microbe_count) %>%
  filter(hap_vap_cap %in% c("HAP", "VAP"))

RA_filt <- fread("results/amr_out/amr_matrices/resfinder.amr.decontam.csv") %>%
  inner_join(meta_filt %>% select(run_id, hap_vap2)) %>%
  select(-hap_vap2) %>%
  column_to_rownames("run_id")

dim(RA_filt)
zeroes <- rownames(RA_filt[rowSums(RA_filt) == 0, ])

# Plot zeroes
mat <- meta_filt %>%
  filter(run_id %in% rownames(RA_filt)) %>%
  mutate(zero = run_id %in% zeroes) %>%
  group_by(hap_vap2) %>%
  summarise(n_absent = sum(zero),
            n_present = sum(!zero)) %>%
  column_to_rownames("hap_vap2")

fisher.test(mat)

RA_filt <- RA_filt[rowSums(RA_filt) > 0, ]

# Bray-Curtis PCoA
bc <- vegdist(as.matrix(RA_filt), method = "jaccard")

PCOA <- pcoa(bc)

eigenvals <- round(PCOA$values$Relative_eig * 100, 1)

plot_df <- as.data.frame(PCOA$vectors) %>%
  rownames_to_column("run_id") %>%
  left_join(meta_filt)

min(plot_df$Axis.1)
max(plot_df$Axis.1)
min(plot_df$Axis.2)
max(plot_df$Axis.2)

p1 <- plot_df %>%
  ggplot(aes(Axis.1, Axis.2, fill = hap_vap2, shape = hap_vap2)) +
  geom_point(color = "black",
             alpha = 0.8, size = 3) +
  scale_fill_manual(values = c("goldenrod", "steelblue", "indianred")) +
  scale_shape_manual(values = c(21, 22, 23)) +
  theme_classic() + 
  ylim(-0.6, 0.7) +
  xlim(-0.7, 0.7) +
  labs(x = str_glue("PCo1 ({eigenvals[1]}%)"),
       y = str_glue("PCo2 ({eigenvals[2]}%)")) +
  theme(legend.position = "none",
        axis.title = element_blank())

a1 <- plot_df %>%
  ggplot(aes(y = Axis.1, x = hap_vap2, fill = hap_vap2)) +
  geom_boxplot() +
  scale_fill_manual(values = c("goldenrod", "steelblue", "indianred")) +
  theme_bw() +
  geom_pwc(label.size = 3, step.increase = 0.05) +
  coord_flip() +
  ylim(-0.7, 0.7) +
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_blank()) +
  labs(y = str_glue("PCo1 ({eigenvals[1]}%)"))

a2 <- plot_df %>%
  ggplot(aes(x = hap_vap2, y = Axis.2, fill = hap_vap2)) +
  geom_boxplot() +
  scale_fill_manual(values = c("goldenrod", "steelblue", "indianred")) +
  theme_bw() +
  geom_pwc(label.size = 3, step.increase = 0.05) +
  ylim(-0.6, 0.7) +
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank()) +
  labs(y = str_glue("PCo2 ({eigenvals[2]}%)"))

ggarrange(a2, p1, ggplot(), a1, nrow = 2, ncol = 2,
          heights = c(2.5, 1), widths = c(1, 4),
          align = "hv")

# ggsave("results/metagenomic_out/PCoA.decontamination.genus.pdf",
#        height = 4, width = 6)


# Logistic regression
# V-HAP versus HAP
# vent_nonvent <- glm(ventilation ~ recent_abx + sample_type + microbial_reads + Axis.1 + Axis.2,
#                     data = plot_df,
#                     family = "binomial")
# 
# summary(vent_nonvent)
# 
# hap_vhap <- glm(ventilation ~ recent_abx + sample_type + microbial_reads + Axis.1 + Axis.2,
#                 data = plot_df %>% filter(hap_vap2 %in% c("V-HAP", "NV-HAP")),
#                 family = "binomial")
# summary(hap_vhap)
# 
# hap_vap <- glm(ventilation ~ recent_abx + sample_type + microbial_reads + Axis.1 + Axis.2 + Axis.3,
#                data = plot_df %>% filter(hap_vap2 %in% c("VAP", "NV-HAP")),
#                family = "binomial")
# 
# summary(hap_vap)
linreg1 <- lm(Axis.1 ~ hospital + sample_type + log10(microbial_reads) + hap_vap2,
             data = plot_df)

anova(linreg1)

linreg2 <- lm(Axis.2 ~ hospital + sample_type + log10(microbial_reads) + hap_vap2,
              data = plot_df)
anova(linreg2)
hist(linreg$residuals)

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
#   left_join(meta_filt) %>%
#   ggplot(aes(x = PC1, y = PC2, color = hap_vap2, shape = hap_vap2)) +
#   geom_point() +
#   # geom_text(aes(label = run_id)) +
#   scale_color_manual(values = c("indianred", "blue", "goldenrod", "grey50", "olivedrab"))+
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