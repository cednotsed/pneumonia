rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)

meta <- fread("data/metadata/parsed_patient_metadata.filt.csv")
parsed <- fread("results/amr_out/amr_matrices/resfinder.amr.decontam.csv") %>%
  filter(run_id %in% meta$run_id)

mat <- parsed %>%
  pivot_longer(!run_id, names_to = "abx", values_to = "presence") %>%
  mutate(abx = ifelse(abx %in% c("pristinamycin_iia", "dalfopristin", "virginiamycin_m"),
                      "streptogramin_a", abx)) %>%
  mutate(abx = ifelse(abx %in% c("pristinamycin_ia", "quinupristin", "virginiamycin_s"),
                      "streptogramin_b", abx)) %>%
  mutate(abx = ifelse(abx %in% c("ampicillin_clavulanic_acid", "amoxicillin_clavulanic_acid", "virginiamycin_m"),
                      "amoxicillin_clavulanic acid", abx)) %>%
  mutate(abx = ifelse(abx %in% c("amoxicillin", "ampicillin"),
                      "aminopenicillin", abx)) %>%
  mutate(abx = ifelse(abx %in% c("clindamycin", "lincomycin"),
                      "clindamycin/lincomycin", abx)) %>%
  mutate(abx = ifelse(abx %in% c("doxycycline", "tetracycline"),
                      "tetracycline", abx)) %>%
  group_by(run_id, abx) %>%
  summarise(n = sum(presence)) %>%
  mutate(presence = ifelse(n > 0, 1, 0)) %>% 
  select(-n) %>%
  pivot_wider(id_cols = run_id, names_from = abx, values_from = presence)

mat %>%
  fwrite("results/amr_out/amr_matrices/resfinder.amr.decontam.agg.csv")
