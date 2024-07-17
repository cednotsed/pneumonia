rm(list = ls())
setwd("C:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)
require(ggpubr)
require(scales)
require(ggsci)
require(ggforce)

meta <- fread("data/metadata/parsed_patient_metadata.filt.csv")

meta_filt <- meta %>%
  filter(hap_vap_cap %in% c("HAP", "VAP", "CAP"))

control_meta <- meta %>%
  filter(hap_vap_cap == "Water control")

micro_meta <- fread("data/metadata/parsed_microbiology_results.csv")

bact_meta <- fread("data/metadata/parsed_microbiology_results.bacterial_sp_only.csv") %>%
  filter(species_resolved)

# Number of invalid tests
micro_meta %>% 
  filter(validity == "Invalid") %>%
  group_by(method) %>%
  summarise(n_failure = n_distinct(run_id),
            perc_failure = n_distinct(run_id) / nrow(micro_meta) * 100)

# No. patients with identifiable pathogens
micro_meta_filt <- micro_meta %>% 
  filter(validity == "Valid",
         bugs != "Negative")

micro_meta_filt %>%
  group_by(run_id) %>%
  summarise(n = n()) %>%
  arrange(n)

# No. of patients with viral/bacterial infections
infection_counts <- micro_meta_filt %>%
  group_by(run_id) %>%
  summarise(n_virus = sum(grepl("virus", bugs, ignore.case = T)),
            n_bact = sum(!grepl("virus", bugs, ignore.case = T))) %>%
  mutate(infection_type = case_when(n_virus > 0 & n_bact == 0 ~ "Viral only",
                                    n_virus > 0 & n_bact > 0 ~ "Bacterial and viral",
                                    n_virus == 0 & n_bact > 0 ~ "Bacterial only",
                                    TRUE ~ NA))

plot_df <- infection_counts %>%  
  group_by(infection_type) %>%
  summarise(n = n()) %>%
  mutate(prop = n / sum(n))

undiagnosed_ids <- deframe(micro_meta %>%
  filter(!(run_id %in% infection_counts$run_id)) %>%
  distinct(run_id))

length(undiagnosed_ids)

# Save undiagnosed pathogens
tibble(x = undiagnosed_ids) %>%
  fwrite("results/benchmarking_out/undiagnosed_ids.txt", col.names = F)

# Plot piechart
plot_df %>% 
  ggplot(aes(x = 1, y = n, fill = infection_type)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = n, hjust = 0)) +
  coord_polar("y") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  scale_fill_manual(limits = c("Bacterial only", "Viral only", "Bacterial and viral"),
                    values = c("indianred", "deepskyblue4", "darkorchid4")) +
  labs(fill = "Infection type")

ggsave("results/benchmarking_out/bacterial_viral.pdf", dpi = 600, width = 5, height = 5)

# Polymicrobial versus monomicrobial
poly_df <- micro_meta_filt %>%
  mutate(bugs = ifelse(grepl("Enterobacter", bugs), "Enterobacter spp.", bugs)) %>%
  mutate(bugs = ifelse(grepl("Acinetobacter", bugs), "Acinetobacter spp.", bugs)) %>%
  mutate(bugs = ifelse(grepl("Proteus", bugs), "Proteus spp.", bugs)) %>%
  group_by(run_id) %>%
  summarise(n = n_distinct(bugs))

poly_df %>%  
  arrange(desc(n)) %>%
  group_by(n) %>%
  summarise(n_patients = n()) %>%
  ggplot(aes(x = n, y = n_patients)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black") +
  labs(x = "No. microbes identified", y = "No. patients") +
  theme_classic()

ggsave("results/benchmarking_out/polymicrobial.pdf", width = 3, height = 3)
poly_df %>%
  mutate(poly = n > 1) %>%
  group_by(poly) %>%
  summarise(n = n()) %>%
  mutate(prop = n / sum(n))
  
  
# sequence_counts <- long_df %>%
#   group_by(run_id) %>%
#   summarise(n = sum(rel_a > 0)) %>%
#   mutate(method = "sequencing")
# 
# merged <- bind_rows(sequence_counts, micro_counts) %>%
#   complete()
# 
# merged %>% View()
# merged %>%
#   ggplot(aes(x = method, y = n)) +
#   geom_boxplot()
# 
# micro_meta %>% 
#   group_by(run_id, method) %>%
#   summarise(n = n()) %>%
#   group_by(run_id)
#   ggplot(aes(x = method, y = n(), fill = factor(n))) +
#   geom_boxplot()
#   left_joi
