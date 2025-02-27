setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(vegan)
require(hillR)
require(ggpubr)

meta <- fread("data/metadata/parsed_patient_metadata.filt.csv") %>%
  filter(hap_vap_cap %in% c("HAP", "VAP"))

RA_filt <- fread("results/tax_classification_out/abundance_matrices/RA.G.zeroed.decontam.2.csv") %>%
  filter(run_id %in% meta$run_id) %>%
  column_to_rownames("run_id")

# Remove zero taxa
RA_filt <- RA_filt[rowSums(RA_filt) != 0, ]

hill_shannon <- hill_taxa(comm = RA_filt, q = 1)

hshan_df <- tibble(run_id = names(hill_shannon),
                   diversity = hill_shannon) %>%
  left_join(meta) %>%
  mutate(recent_abx = gsub("sample", "sampling", recent_abx)) %>%
  mutate(recent_abx = factor(recent_abx, c(NA, "After sampling", "Before sampling")))

count_df <- hshan_df %>%
  group_by(hap_vap2) %>%
  summarise(n = n())

hshan_df %>%
  ggplot(aes(x = hap_vap2, y = diversity, fill = hap_vap2)) +
  geom_boxplot() +
  geom_pwc() +
  geom_text(aes(x = hap_vap2, y = 0, label = str_glue("n={n}")),
            data = count_df,
            position = position_dodge(width = 0.8)) +
  labs(x = "Patient group", 
       y = "Hill-Shannon diversity",
       fill = "Pneumonia subtype") +
  scale_fill_manual(values = c("goldenrod", "steelblue", "indianred")) +
  theme_classic() +
  theme(legend.position = "none")

ggsave("results/metagenomic_out/microbial_diversity.genus.pdf", dpi = 600, 
       width = 5, height = 3)

# # V-HAP versus HAP
# vent_nonvent <- glm(ventilation ~ recent_abx + sample_type + microbial_reads + diversity,
#                   data = hshan_df,
#                   family = "binomial")
# summary(vent_nonvent)
# 
# hap_vhap <- glm(ventilation ~ recent_abx + sample_type + microbial_reads + diversity,
#                     data = hshan_df %>% filter(hap_vap2 %in% c("V-HAP", "NV-HAP")),
#                     family = "binomial")
# summary(hap_vhap)
# 
# hap_vap <- glm(ventilation ~ recent_abx + sample_type + microbial_reads + diversity,
#                 data = hshan_df %>% filter(hap_vap2 %in% c("VAP", "NV-HAP")),
#                 family = "binomial")
# 
# summary(hap_vap)

linreg <- lm(log10(diversity) ~ hospital + log10(microbial_reads) + hap_vap2,
   data = hshan_df)

anova(linreg)
summary(linreg)
hist(linreg$residuals)
