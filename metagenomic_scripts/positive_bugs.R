setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(doParallel)

RA_filt <- fread("results/metagenomic_out/RA.S.zeroed.csv")
meta <- fread("data/metadata/parsed_patient_metadata.csv")

# Get bug list
bug_list <- unique(c(meta$micro_organism, meta$biofire_organism, meta$curetis_organism))
bug_list <- gsub("\\s*\\([^\\)]+\\)","", bug_list)
bug_list <- unique(unlist(str_split(bug_list, ";")))
bug_list <- bug_list[!(bug_list %in% c("Negative", "", "Invalid test"))]
bug_list <- bug_list[!grepl("complex|group|Coliform|spp|virus|-|Pneumocystis", bug_list)]

morsels <- foreach(bug_name = bug_list) %do% {
  # bug_name <- bug_list[1]
  positive_patients <- meta %>%
    filter(grepl(bug_name, micro_organism)|
             grepl(bug_name, biofire_organism)|
             grepl(bug_name, curetis_organism))
  
  water <- meta %>%
    filter(hap_vap_cap == "Water control")
  
  taxa_list <- colnames(RA_filt)[grepl(bug_name, colnames(RA_filt))]
  
  RA_filt %>%
    select(all_of(c("run_id", bug_name))) %>%
    mutate(sample_type = case_when(run_id %in% positive_patients$run_id ~ "Positive",
                                   run_id %in% water$run_id ~ "Blank",
                                   TRUE ~ "Negative")) %>%
    mutate(sample_type = factor(sample_type, c("Negative", "Positive", "Blank"))) %>%
    pivot_longer(!c("run_id", "sample_type"), names_to = "bug", values_to = "rel_a")
}

bind_rows(morsels) %>%
  ggplot(aes(x = sample_type, y = rel_a, fill = bug)) +
  geom_boxplot() +
  facet_wrap(.~bug) +
  labs(x = "Sample type", y = "Relative abundance") +
  theme_bw() +
  theme(legend.position = "none")

