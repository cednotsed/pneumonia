setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(vegan)
require(hillR)
require(ggpubr)
require(randomcoloR)

meta <- fread("data/metadata/parsed_patient_metadata.filt.csv")

file_dir <- "results/amr_out/resfinder_out/all_phenotypes/"
file_list <- list.files(file_dir, full.names = T)

id_list <- foreach(file = file_list, .combine = "c") %do% {
  id <- gsub(file_dir, "", file)
  
  return(id)
}

res_df <- foreach(file = file_list, .combine = "bind_rows") %do% {
  id <- gsub(file_dir, "", file)
  
  fread(file) %>%
    mutate(run_name = id)
}

parsed <- res_df %>%
  rename_all(~tolower(gsub(" |-", "_", .x))) %>%
  rename(antimicrobial = `#_antimicrobial`) %>%
  mutate(run_id = gsub("INHALE_FRESH_", "", run_name)) %>%
  mutate(run_id = gsub(".pheno_table.txt", "", run_id)) %>%
  mutate(run_id = gsub("barcode", "", run_id)) %>%
  mutate(run_id = gsub("-0", "-", run_id)) %>%
  mutate(run_id = gsub("-", "_", run_id)) %>%
  mutate(run_id = gsub("12a", "12", run_id)) %>%
  inner_join(meta %>% select(run_id, hap_vap_cap)) %>%
  arrange(desc(antimicrobial))

drug_list <- deframe(parsed %>%
                       distinct(antimicrobial) %>%
                       arrange(desc(antimicrobial)))
parsed %>%
  mutate(antimicrobial = factor(antimicrobial, drug_list)) %>%
  ggplot(aes(x = run_id, y = antimicrobial, fill = wgs_predicted_phenotype)) +
  geom_tile() +
  scale_fill_manual(values = c("steelblue", "red")) +
  facet_grid(cols = vars(hap_vap_cap), scales = "free", space = "free") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(x = "Sample", y = "Antimicrobial",
       fill = "Predicted phenotype")

ggsave("results/amr_out/amr_profiles.pdf", width = 15, height = 15)

parsed %>%
  group_by(run_id) %>%
  filter(wgs_predicted_phenotype == "Resistant") %>%
  summarise(n = n_distinct(class)) %>%
  ggplot(aes(x = n)) +
  geom_histogram() +
  labs(x = "No. antimicrobial resistance classes", 
       y = "No. patients")

plot_df <- parsed %>%
  group_by(antimicrobial) %>%
  summarise(prop = sum(wgs_predicted_phenotype == "Resistant") / n()) %>%
  arrange(desc(prop)) %>% 
  head(20) %>%
  mutate(antimicrobial = Hmisc::capitalize(antimicrobial))
pal <- distinctColorPalette(n_distinct(plot_df$antimicrobial))
plot_df %>%
  filter(antimicrobial %in% plot_df$antimicrobial) %>%
  mutate(antimicrobial = factor(antimicrobial, plot_df$antimicrobial),
         prop = round(prop, 2)) %>%
  ggplot(aes(x = antimicrobial, y = prop, fill = antimicrobial)) +
  geom_bar(stat = "identity",
           color = "black") +
  # geom_text(aes(label = prop)) +
  scale_fill_manual(values = pal) + 
  labs(x = "Antimicrobial", y = "Prop. of patients resistant to") +
  theme_classic() + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("results/amr_out/amr_barplots.pdf", dpi = 600, width = 5, height = 5)  
ggsave("results/amr_out/amr_barplots.png", dpi = 600, width = 5, height = 5)  
filter(wgs_predicted_phenotype == "Resistant") %>%
  distinct(run_id)
  
  summarise(n = n_distinct(antimicrobial))
id_parsed <- gsub("INHALE_FRESH_|.pheno_table.txt|barcode", "", id_list)
id_parsed <- gsub("-0", "-", id_parsed)
id_parsed <- gsub("-", "_", id_parsed)
id_parsed <- gsub("12a", "12", id_parsed)

id_parsed
  left_join(meta %>% select(run_id, hap_vap_cap)) 
meta %>%
  filter(run_id %in% id_parsed)
parsed %>%
  distinct()
