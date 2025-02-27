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

meta <- fread("data/metadata/parsed_patient_metadata.filt.csv") %>%
  filter(hap_vap_cap %in% c("HAP", "VAP"))

amr_df <- fread("results/amr_out/resfinder_out/resfinder.assemblies.parsed.csv") %>% 
  filter(wgs_predicted_phenotype == "Resistant") %>%
  filter(run_id %in% meta$run_id) %>%
  mutate(antimicrobial = capitalize(antimicrobial)) %>%
  mutate(class = capitalize(class))

morsels <- foreach(idx = seq(nrow(amr_df))) %do% {
  row <- amr_df[idx, ]
  genes <- str_split(row$genetic_background, "\\,\\ ")[[1]]
  genes <- str_split(genes, "\\ ", simplify = T)[, 1]
  
  tibble(genome_name = row$genome_name, run_id = row$run_id, antimicrobial = row$antimicrobial, class = row$class, gene = genes)
}

merged <- bind_rows(morsels) %>%
  left_join(amr_df %>% distinct(genome_name, family, genus, 
                                species, checkm2_completeness, checkm2_contamination, 
                                genome_size, gtdbtk_warnings)) %>%
  mutate(parsed_gene = case_when(grepl("aac", gene) ~ "aac",
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
                                grepl("cepA", gene) ~ "cepA",
                                TRUE ~ gene))

merged %>%
  fwrite("results/amr_out/amr_matrices/resfinder.assemblies.gene_content.csv")

