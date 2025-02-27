rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(ggpubr)
require(viridis)
require(Hmisc)

meta <- fread("data/metadata/parsed_patient_metadata.filt.csv") %>%
  filter(high_microbe_count) %>%
  filter(hap_vap_cap %in% c("HAP", "VAP"))

df <- fread("results/amr_out/amr_matrices/resfinder.amr.decontam.agg.csv") %>%
  pivot_longer(!run_id, names_to = "antimicrobial", values_to = "presence") %>%
  inner_join(meta %>% select(run_id, hap_vap_cap, hap_vap2)) %>%
  mutate(antimicrobial = capitalize(antimicrobial))

class_meta <- fread("results/amr_out/resfinder_out/resfinder.parsed.decontam.csv") %>%
  distinct(antimicrobial, class) %>%
  mutate(antimicrobial = gsub("\\+|\\-", "_", antimicrobial)) %>%
  mutate(antimicrobial = gsub("\\ ", "_", antimicrobial)) %>%
  mutate(antimicrobial = capitalize(antimicrobial)) %>%
  mutate(class = capitalize(class))

count_df <- df %>%
  group_by(hap_vap2) %>%
  summarise(n_total = n_distinct(run_id))

to_keep <- df %>%
  group_by(antimicrobial) %>%
  summarise(n = sum(presence == 1)) %>%
  arrange(desc(n)) %>%
  filter(n >= 10)

plot_df <- df %>%
  group_by(hap_vap2, antimicrobial) %>%
  summarise(n = sum(presence == 1)) %>%
  arrange(desc(n)) %>%
  left_join(count_df) %>%
  mutate(prop = n / n_total) %>%
  left_join(class_meta) %>%
  filter(antimicrobial %in% to_keep$antimicrobial) %>%
  mutate(class = case_when(antimicrobial == "Streptogramin_b" ~ "Streptogramin",
                           antimicrobial == "Streptogramin_a" ~ "Streptogramin",
                           antimicrobial == "Clindamycin/lincomycin" ~ "Lincosamide",
                           antimicrobial == "Aminopenicillin" ~ "Beta-lactam",
                           antimicrobial == "Amoxicillin_clavulanic acid" ~ "Beta-lactam",
                           TRUE ~ class)) %>%
  filter(!(antimicrobial %in% c("Dibekacin", "Fortimicin", "Isepamicin", 
                                "Butirosan", "Lividomycin", "Paromycin", 
                                "Ribostamicin", "Sisomycin"))) %>%
  filter(antimicrobial != "Unknown_beta_lactam")

order_df <- plot_df %>%
  filter(hap_vap2 == "NV-HAP") %>%
  left_join(class_meta) %>%
  arrange(desc(prop))

plot_df %>% 
  mutate(antimicrobial = factor(antimicrobial, unique(order_df$antimicrobial))) %>%
  mutate(class = factor(class, unique(order_df$class))) %>%
  ggplot(aes(x = antimicrobial, y = hap_vap2, fill = prop)) +
  geom_tile(color = "white") +
  facet_grid(cols = vars(class), space = "free", scales = "free") +
  scale_fill_viridis_c() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank()) +
  labs(x = "Predicted resistance", y = "HAP subtype", fill = "Prop. patients") 

ggsave("results/amr_out/amr_profiles.pdf", width = 12, height = 4, dpi = 600)
