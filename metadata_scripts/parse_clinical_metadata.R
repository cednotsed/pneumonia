rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)

meta <- fread("data/metadata/parsed_patient_metadata.filt.csv")

clin_meta <- fread("data/metadata/clinical_metadata.csv")

merged_df <- clin_meta %>%
  left_join(meta) %>%
  filter(!is.na(n_reads)) %>%
  mutate(parsed_outcome = case_when(day_21_outcome %in% c("No data", "Not Known",
                                                          "Patient died of unrelated causes") ~ "Unknown",
                                    day_21_outcome == "Pneumonia cured" ~ "Cured",
                                    day_21_outcome %in% c("Pneumonia not cured") ~ "Not cured",
                                    day_21_outcome == "Patient died of pneumonia" ~ "Died"))

table(merged_df$parsed_outcome)
merged_df %>% 
  group_by(day_21_outcome) %>%
  summarise(n = n())

merged_df %>%
  fwrite("data/metadata/parsed_clinical_metadata.csv")

