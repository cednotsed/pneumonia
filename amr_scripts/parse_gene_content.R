rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(ggpubr)
require(viridis)
require(scales)
require(ggsci)
require(ggforce)
require(Hmisc)

cor_df <- fread("../ARG_host_jumps/results/gene_classification_out/gene_cluster_meta.csv")
cor_df %>% View()
meta <- fread("data/metadata/parsed_patient_metadata.filt.csv") %>%
  filter(hap_vap_cap %in% c("HAP", "VAP"))

amr_df <- fread("results/amr_out/resfinder_out/resfinder.parsed.decontam.csv") %>% 
  filter(wgs_predicted_phenotype == "Resistant") %>%
  filter(run_id %in% meta$run_id) %>%
  mutate(antimicrobial = capitalize(antimicrobial)) %>%
  mutate(class = capitalize(class))

morsels <- foreach(idx = seq(nrow(amr_df))) %do% {
  row <- amr_df[idx, ]
  genes <- str_split(row$genetic_background, "\\,\\ ")[[1]]
  genes <- str_split(genes, "\\ ", simplify = T)[, 1]
  
  tibble(run_id = row$run_id, antimicrobial = row$antimicrobial, class = row$class, gene = genes)
}

merged <- bind_rows(morsels) %>%
  mutate(gene = case_when(grepl("aac", gene) ~ "aac",
                          grepl("blaOXA", gene) ~ "blaOXA",
                          grepl("aad", gene) ~ "aad",
                          grepl("tet", gene) ~ "tet",
                          grepl("qnr", gene) ~ "qnr",
                          grepl("erm", gene) ~ "erm",
                          grepl("blaTEM", gene) ~ "blaTEM",
                          grepl("blaSHV", gene) ~ "blaSHV",
                          grepl("blaOXY", gene) ~ "blaOXY",
                          grepl("blaCTX", gene) ~ "blaCTX",
                          grepl("blaBRO", gene) ~ "blaBRO",
                          grepl("blaCMY", gene) ~ "blaCMY",
                          grepl("blaDHA", gene) ~ "blaDHA",
                          grepl("vga", gene) ~ "vga",
                          grepl("ant", gene) ~ "ant",
                          grepl("aph", gene) ~ "aph",
                          grepl("sul", gene) ~ "sul",
                          grepl("Oqx", gene) ~ "Oqx",
                          grepl("msr", gene) ~ "msr",
                          grepl("Nar", gene) ~ "Nar",
                          grepl("cat", gene) ~ "cat",
                          grepl("mph", gene) ~ "mph",
                          grepl("fos", gene) ~ "fos",
                          grepl("lsa", gene) ~ "lsa",
                          grepl("fus", gene) ~ "fus",
                          grepl("dfr", gene) ~ "dfr",
                          grepl("cfx", gene) ~ "cfx",
                          TRUE ~ gene))

merged %>%
  fwrite("results/amr_out/amr_matrices/resfinder.gene_content.csv")

