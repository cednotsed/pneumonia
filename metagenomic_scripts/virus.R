setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(vegan)
require(hillR)
require(ggpubr)

meta <- fread("data/metadata/parsed_patient_metadata.csv")

virus_patients <- meta %>%
  filter(grepl("virus", curetis_organism, ignore.case = T) |
           grepl("virus", biofire_organism, ignore.case = T))

virus_patients %>% View()

tax_meta <- fread("data/metadata/k2_database_taxonomy.csv")
tax_meta$group
RA_filt <- fread("results/metagenomic_out/abundance_matrix.S.tsv")
RA_filt[1:5, 1:5]
dat %>%
  rownames_to_column("run_id") %>%
  pivot_longer(!run_id, 
               names_to = "taxa", 
               values_to = "rel_a") %>%
  filter(rel_a != 0) %>%
  filter(grepl("virus|influ", taxa)) %>% View()

# Get bug list
bug_list <- unique(c(meta$micro_organism, meta$biofire_organism, meta$curetis_organism))
bug_list <- gsub("\\s*\\([^\\)]+\\)","", bug_list)
bug_list <- unique(unlist(str_split(bug_list, ";")))
bug_list <- bug_list[!(bug_list %in% c("Negative", "", "Invalid test"))]
bug_list <- bug_list[!grepl("complex|group|Coliform|spp|virus|-|Pneumocystis", bug_list)]

RA_filt <- fread("results/metagenomic_out/RA.S.zeroed.csv") %>%
  column_to_rownames("run_id")

plot_df <- RA_filt %>%
  rownames_to_column("run_id") %>%
  pivot_longer(!run_id, names_to = "taxa", values_to = "rel_a") %>%
  filter(rel_a != 0) %>%
  mutate(is_pathogen = ifelse(taxa %in% bug_list, "pathogen", "nonpathogen")) %>%
  group_by(run_id, is_pathogen) %>%
  summarise(sum_a = sum(rel_a)) %>%
  pivot_wider(id_cols = run_id, names_from = "is_pathogen", values_from = "sum_a") %>%
  mutate(pathogen = replace_na(pathogen, 0),
         nonpathogen = replace_na(nonpathogen, 0)) %>%
  mutate(path_ratio = pathogen / nonpathogen) %>%
  left_join(meta) %>%
  filter(!(hosp_los_hours %in% c("", "Not Known"))) %>%
  mutate(hosp_los_hours = as.numeric(hosp_los_hours)) %>%
  filter(pathogen != 0 & nonpathogen != 0)

plot_df %>%
  filter(pathogen != 0 & nonpathogen != 0) %>%
  ggplot(aes(x = hosp_los_hours, y = path_ratio)) +
  geom_point() +
  geom_smooth()

g1 <- glm(path_ratio ~ sample_type + post_pcr_qubit + hosp_los_hours,
          data = plot_df,
          family = Gamma(link = "inverse"))

summary(g1)

plot_df %>%
  filter(hap_vap_cap == "HAP") %>%
  filter(hosp_los_hours < 100) %>%
  ggplot(aes(x = as.numeric(hosp_los_hours), y = pathogen)) +
  geom_point()

meta %>%
  
  ggplot(aes(x = log10(as.numeric(hosp_los_hours)))) +
  geom_histogram()
ggplot(aes(x = vent_cat, y = max_a)) +
  geom_boxplot()
geom_point() +
  geom_smooth(method = "lm")
# Remove zero taxa
RA_filt <- RA_filt[rowSums(RA_filt) != 0, ]

hill_shannon <- hill_taxa(comm = RA_filt, q = 1)

hshan_df <- tibble(run_id = names(hill_shannon),
                   diversity = hill_shannon) %>%
  left_join(meta) %>%
  mutate(recent_abx = gsub("sample", "sampling", recent_abx)) %>%
  mutate(recent_abx = factor(recent_abx, c(NA, "After sampling", "Before sampling"))) %>%
  filter(!(hosp_los_hours %in% c("", "Not Known")))

hshan_df %>%
  filter(hap_vap_cap == "VAP") %>%
  ggplot(aes(x = log10(vent_lenght_hours), y = diversity, color = post_pcr_qubit)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "Time on ventilator (hrs)", y = "Hill-Shannon diversity")

meta_filt$vent_lenght_hours

lr <- lm(diversity ~ sample_type + log10(post_pcr_qubit) + vent_lenght_hours,
         data = hshan_df %>%
           filter(hap_vap_cap == "VAP") %>%
           filter(post_pcr_qubit > 10))
summary(lr)

hshan_df %>%
  ggplot(aes(x = log10(post_pcr_qubit), y = diversity)) +
  geom_point()
View(meta)

hshan_df