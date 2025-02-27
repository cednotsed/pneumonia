setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(vegan)
require(hillR)
require(ggpubr)
require(ghibli)

meta <- fread("data/metadata/parsed_patient_metadata.filt.csv") %>%
  left_join(fread("data/metadata/antibiotic_metadata.csv"))

meta_filt <- meta %>%
  filter(!is.na(antibiotics))

controls <- meta %>%
  filter(hap_vap_cap == "Healthy")

# fread("data/metadata/parsed_clinical_metadata.csv") %>% View()
tax_meta <- fread("databases/k2_pluspfp_20240112/inspect.txt")
fungi <- tax_meta[501:801, ]$V6
bacteria <- tax_meta[996:20619, ]$V6

# Remove fungi
read_filt <- fread("results/tax_classification_out/abundance_matrices/read_counts.G.zeroed.decontam.2.csv") %>%
  filter(run_id %in% c(meta_filt$run_id, controls$run_id)) %>%
  column_to_rownames("run_id") %>%
  select(any_of(bacteria))

# Rescale relative abundance
otu_to_RA <- function(df) {
  mat <- as.matrix(df)
  RA_df <- as.data.frame(mat / rowSums(mat))
  RA_df[is.na(RA_df)] <- 0
  colnames(RA_df) <- colnames(df)
  
  return(RA_df)
}

bact_RA <- otu_to_RA(read_filt)
bact_RA <- bact_RA[ , colSums(bact_RA) != 0]

# Remove zero taxa
hill_shannon <- hill_taxa(comm = bact_RA, q = 1)

hshan_df <- tibble(run_id = names(hill_shannon),
                   diversity = hill_shannon) %>%
  inner_join(meta) %>%
  mutate(antibiotics = ifelse(is.na(antibiotics), "Healthy", antibiotics)) %>%
  mutate(antibiotics = factor(antibiotics, c("Healthy", "No antibiotics", 
                                             "Narrow spectrum", "Broad spectrum")))

count_df <- hshan_df %>%
  group_by(antibiotics) %>%
  summarise(n = n())

pal <- ghibli_palettes$KikiLight[c(3, 4, 6, 7)]

hshan_df %>%
  ggplot(aes(x = antibiotics, y = diversity, fill = antibiotics)) +
  geom_boxplot() +
  geom_pwc() +
  geom_text(aes(x = antibiotics, y = 0, label = str_glue("n={n}")),
            data = count_df,
            position = position_dodge(width = 0.8)) +
  labs(x = "Patient group", 
       y = "Hill-Shannon Diversity (Richness)") +
  scale_fill_manual(values = pal) +
  theme_classic() +
  theme(legend.position = "none")
        # axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("results/clinical_out/diversity_vs_abx.pdf", dpi = 600, height = 3, width = 5)


hshan_df %>%
  mutate(group = ifelse(antibiotics == "healthy", "healthy", "patients")) %>%
  ggplot(aes(x = group, y = diversity)) +
  geom_boxplot() +
  geom_pwc()
  
# hshan_df %>%
#   filter(is.na(antibiotics))
# 
# 
# hshan_df %>%
#   ggplot(aes(x = hap_vap_cap, y = diversity)) +
#   geom_boxplot() +
#   geom_pwc() +
#   geom_text(aes(x = hap_vap_cap, y = 0, label = str_glue("n={n}")),
#             data = count_df,
#             position = position_dodge(width = 0.8)) +
#   labs(x = "Patient group", 
#        y = "Hill-Shannon Diversity (Richness)",
#        fill = "Antibiotics given") +
#   theme_bw() +
#   scale_fill_discrete(na.value = "grey")
# 
# 
# 
# shannon_lr <- glm((as.numeric(factor(hap_vap_cap)) - 1) ~ recent_abx + diversity,
#                   data = hshan_df,
#                   family = "binomial")
# 
# summary(shannon_lr)
# 
# hill_simpson <- hill_taxa(comm = RA_filt, q = 2)
# 
# hsimp_df <- tibble(run_id = names(hill_simpson),
#                    diversity = hill_simpson) %>%
#   left_join(meta)
# # filter(hap_vap_cap != "CAP")
# 
# count_df <- hsimp_df %>%
#   group_by(hap_vap_cap, recent_abx) %>%
#   summarise(n = n())
# 
# hsimp_df %>% 
#   # mutate(hap_vap_cap = ifelse(hap_vap_cap == "Water control", hap_vap_cap, "Pneumonia")) %>%
#   ggplot(aes(x = hap_vap_cap, y = log10(diversity))) +
#   geom_violin() +
#   geom_pwc()
# geom_text(aes(x = hap_vap_cap, y = 1, label = str_glue("n={n}")),
#           data = count_df,
#           position = position_dodge(width = 0.8)) +
#   labs(x = "Patient group", 
#        y = "Hill-Simpson Diversity (Evenness)",
#        fill = "Antibiotics given")
# 
# hsimp_df %>% 
#   filter(is.na(hap_vap_cap))
# shannon_lr <- glm((as.numeric(factor(hap_vap_cap)) - 1) ~ recent_abx + diversity,
#                   data = hshan_df,
#                   family = "binomial")
# 
# count_df <- hsimp_df %>%
#   group_by(hap_vap_cap, recent_abx, sample_type) %>%
#   summarise(n_samples = n()) %>%
#   ungroup() %>%
#   complete(hap_vap_cap, recent_abx, sample_type) %>%
#   mutate(n_samples = replace_na(n_samples, 0))
# 
# # hsimp_df %>%
# #   ggplot(aes(x = hap_vap_cap, y = diversity, fill = recent_abx)) +
# #   geom_boxplot() +
# #   facet_grid(.~sample_type) +
# #   geom_text(aes(x = hap_vap_cap, 
# #                 y = Inf, 
# #                 label = str_glue("n={n_samples}"),
# #                 color = recent_abx),
# #             position = position_dodge(width = 0.9),
# #             vjust = 2,
# #             size = 2,
# #             data = count_df) +
# #   geom_pwc()
# 
# hsimp_df %>%
#   filter(hap_vap_cap %in% c("Water control", "VAP")) %>%
#   # mutate(hap_vap_cap = ifelse(hap_vap_cap == "Water control", hap_vap_cap, "Pneumonia")) %>%
#   group_by(hap_vap_cap) %>%
#   summarise(mean_div = mean(diversity),
#             sd = sd(diversity),
#             lower_q = quantile(diversity, 0.25),
#             upper_q = quantile(diversity, 0.75),
#             min = min(diversity),
#             max = max(diversity))
# 
# test <- hsimp_df %>%
#   filter(hap_vap_cap == "Water control")
# 
# shapiro.test(test$diversity)
