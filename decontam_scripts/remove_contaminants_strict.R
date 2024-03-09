setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(randomcoloR)

RA_filt <- fread("results/metagenomic_out/RA.G.zeroed.csv")
meta <- fread("data/metadata/parsed_patient_metadata.csv")

meta <- fread("data/metadata/parsed_patient_metadata.csv")

# Get bug list
bug_list <- unique(c(meta$micro_organism, meta$biofire_organism, meta$curetis_organism))
bug_list <- gsub("\\s*\\([^\\)]+\\)","", bug_list)
bug_list <- unique(unlist(str_split(bug_list, ";")))
bug_list <- bug_list[!(bug_list %in% c("Negative", "", "Invalid test"))]
bug_list <- bug_list[!grepl("virus|-|Pneumocystis", bug_list)]
bug_list <- unique(str_split(bug_list, "\\ ", simplify = T)[, 1])

# Find contaminant taxa
control_long <- RA_filt %>%
  left_join(meta %>% select(run_id, hap_vap_cap)) %>%
  filter(hap_vap_cap == "Water control") %>%
  select(-hap_vap_cap) %>%
  pivot_longer(!run_id, names_to = "taxa", values_to = "rel_a") %>%
  filter(rel_a != 0) %>%
  group_by(taxa) %>%
  summarise(n = n_distinct(run_id)) %>%
  arrange(desc(n)) %>%
  mutate(is_pathogen = taxa %in% bug_list)

# pal <- distinctColorPalette(n_distinct(control_long$taxa))
control_long %>%
  mutate(taxa = factor(taxa, unique(control_long$taxa))) %>%
  mutate(is_pathogen = ifelse(is_pathogen, "Yes", "No")) %>%
  ggplot(aes(x = taxa, y = n, fill = is_pathogen)) +
  geom_bar(stat = "identity") +
  labs(x = "Contaminant genera",
       y = "Number of blank samples",
       fill = "Identified via culture/PCR?") +
  scale_fill_manual(values = c("#64A6E3", "#E68FAA")) +
  theme_bw() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 90, 
                                   hjust = 1,
                                   vjust = 0.3))

decontam <- RA_filt %>%
  select(-all_of(control_long$taxa))

control_long %>%
  select(taxa) %>% 
  nrow()

fwrite(decontam,
       "results/metagenomic_out/RA.S.zeroed.decontam.csv")

# control_long
# control_filt
# foreach(run_name = unique(meta$run)) %do% {
#   run_name = unique(meta$run)[1]
#   patient_ids <- deframe(meta %>% 
#                            filter(hap_vap_cap != "Water control") %>%                          
#                            filter(run == run_name) %>%
#                            select(run_id))
#   
#   control_ids <- deframe(meta %>% 
#                            filter(hap_vap_cap == "Water control") %>%                          
#                            filter(run == run_name) %>%
#                            select(run_id))
#   
#   RA_filt %>%
#    filter(run_id %in% run_ids)
# }
# 
# 
# control_ids <- meta %>% 
#   filter(hap_vap_cap == "Water control")
# 
# RA_filt %>%
#   left_join(meta %>% select(run_id, run, hap_vap_cap)) %>%
#   filter(hap_vap_cap == "Water control") %>%
#   pivot_longer(!c("hap_vap_cap", "run_id", "run"),
#                names_to = "taxa",
#                values_to = "rel_a") %>%
#   group_by(run_id) %>%
#   filter(rel_a != 0) %>%
#   arrange(desc(rel_a)) %>% View()
# 
# RA_filt %>%
#   filter(run_id == "10_1") %>%
#   pivot_longer(!c("run_id"),
#                names_to = "taxa",
#                values_to = "rel_a") %>%
#   arrange(desc(rel_a)) %>%
#   View()
