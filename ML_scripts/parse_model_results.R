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

file_dir <- "results/ML_out/raw_results/"

file_list <- list.files(file_dir, full.names = T)

morsels <- foreach(file_name = file_list) %do% {
  fread(file_name) %>%
    mutate(id = file_name)
}

df <- bind_rows(morsels) %>%
  mutate(id = gsub(file_dir, "", id)) %>%
  mutate(id = gsub(".csv", "", id)) %>%
  separate(id, c("exp", "type"), "\\.")

df %>%
  ggplot(aes(x = exp, y = test_F1, fill = type)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_manual(values = c("olivedrab")) +
  labs(x = "XGBoost Model", y = "Performance (F1 score)") +
  theme_bw()

ggsave("results/ML_out/model_results.pdf", dpi = 600, height = 4, width = 4)
