rm(list = ls())
setwd("C:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(ggpubr)
require(viridis)

meta <- fread("data/metadata/parsed_patient_metadata.filt.csv")

meta_filt <- meta %>%
  filter(hap_vap_cap %in% c("HAP", "VAP")) %>%
  filter(high_microbe_count) %>%
  mutate(ventilation = ifelse(ventilation, "VAP", "HAP"))

df_filt <- fread("results/tax_classification_out/abundance_matrices/RA.G.zeroed.decontam.2.csv") %>%
  filter(run_id %in% meta_filt$run_id) %>%
  column_to_rownames("run_id")

df_filt[df_filt > 0] <- 1

df_filt <- df_filt %>%
  rownames_to_column("run_id")

merged <- df_filt %>%
  left_join(meta_filt)

logreg <- glm(ventilation ~ recent_abx + microbial_reads + sample_type + Lacticaseibacillus + Klebsiella + 	
                Veillonella + Actinomyces + Pseudomonas + Rothia,
    data = merged %>% 
      mutate(ventilation = ifelse(ventilation == "VAP", 1, 0)),
    family = "binomial")

summary(logreg)
coefficients(summary(logreg))[, "Pr(>|z|)"]
aod <- anova(logreg)
aod["Deviance"] / aod["NULL", "Resid. Dev"] * 100

sputum_logreg <- glm(ventilation ~ recent_abx + microbial_reads + Lacticaseibacillus + Klebsiella + 	
                     Veillonella + Actinomyces + Pseudomonas + Rothia,
                     data = merged %>% 
                       mutate(ventilation = ifelse(ventilation == "VAP", 1, 0)) %>%
                       filter(sample_type == "SPU"),
                     family = "binomial")

summary(logreg)

lac <- glm(Lacticaseibacillus ~ recent_abx + microbial_reads + sample_type + ventilation,
    data = merged %>% 
      mutate(ventilation = ifelse(ventilation == "VAP", 1, 0)),
    family = "binomial")

summary(lac)

kleb <- glm(Klebsiella ~ recent_abx + microbial_reads + sample_type + ventilation,
           data = merged %>% 
             mutate(ventilation = ifelse(ventilation == "VAP", 1, 0)),
           family = "binomial")
summary(kleb)

veil <- glm(Veillonella ~ recent_abx + microbial_reads + sample_type + ventilation,
            data = merged %>% 
              mutate(ventilation = ifelse(ventilation == "VAP", 1, 0)),
            family = "binomial")

summary(veil)

act <- glm(Actinomyces ~ recent_abx + microbial_reads + sample_type + ventilation,
           data = merged %>% 
             mutate(ventilation = ifelse(ventilation == "VAP", 1, 0)),
           family = "binomial")

summary(act)

psu <- glm(Pseudomonas ~ recent_abx + microbial_reads + sample_type + ventilation,
           data = merged %>% 
             mutate(ventilation = ifelse(ventilation == "VAP", 1, 0)),
           family = "binomial")

summary(psu)

roth <- glm(Rothia ~ recent_abx + microbial_reads + sample_type + ventilation,
           data = merged %>% 
             mutate(ventilation = ifelse(ventilation == "VAP", 1, 0)),
           family = "binomial")

summary(roth)
