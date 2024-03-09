rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)

tax_meta <- fread("databases/k2_pluspf_20231009/inspect.txt", skip = 7) %>%
  select(taxid = V5, taxon_name = V6)

meta <- fread("data/metadata/parsed_patient_metadata.filt.csv")
cov_df <- fread("results/irep_out/coverage_results.parsed.tsv")
irep_df <- fread("results/irep_out/irep_results.parsed.tsv")

merged <- cov_df %>%
  left_join(irep_df) %>%
  left_join(tax_meta) %>%
  inner_join(meta %>% select(-n_reads)) %>%
  separate(taxon_name, into = c("genus"), "\\ ")


test <- merged %>%
  # filter(!is.na(bPTR)) %>%
  # group_by(run_id, genus, hap_vap_cap) %>%
  # summarise(max_bPTR = max(bPTR, na.rm = T)) %>%
  filter(genus == "Klebsiella",
         !is.na(bPTR)) %>%
  distinct(run_id)
  ggplot(aes(x = genus, y = run_id, fill = max_bPTR)) +
  geom_tile() +
  facet_grid(rows = vars(hap_vap_cap), scales = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3),
        axis.text.y = element_text(size = 5))

merged %>%
  group_by(hap_vap_cap) %>%
  summarise(n = n_distinct(run_id))
merged %>%
  filter(hap_vap_cap == "Water control") %>% View()

meta %>%
  filter(run_id %in% test$run_id) %>%
  select(run_id, micro_organism, curetis_organism, biofire_organism) %>% View()

merged_filt <- merged %>%
  filter(!is.na(bPTR))

# Get culture genera
culture_res <- deframe(meta %>%
  distinct(micro_organism))

bug_list <- foreach(entry = culture_res, .combine = "c") %do% {
  str_split(entry, "\\;", simplify = T)[1, ]
}

bug_list_filt <- unique(bug_list)
bug_list_filt <- bug_list_filt[bug_list_filt != ""]
bug_list_filt
foreach()

bPTR
# meta %>% 
#   mutate(pre_wash_qubit = gsub("<0.05", "0.0001", pre_wash_qubit)) %>%
#   mutate(pre_wash_qubit = as.numeric(pre_wash_qubit)) %>%
#   ggplot(aes(x = log10(post_pcr_qubit), fill = hap_vap_cap)) +
#   geom_histogram() +
#   facet_grid(rows = vars(hap_vap_cap))
