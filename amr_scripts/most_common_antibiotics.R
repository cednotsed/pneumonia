rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(ggpubr)
require(viridis)
require(Hmisc)

meta <- fread("data/metadata/parsed_patient_metadata.filt.csv")

drug_meta <- fread("results/amr_out/resfinder_out/resfinder.parsed.decontam.csv") %>%
  distinct(antimicrobial, class) %>%
  mutate(antimicrobial = gsub("\\+|\\-", "_", antimicrobial)) %>%
  mutate(antimicrobial = gsub("\\ ", "_", antimicrobial)) %>%
  mutate(antimicrobial = capitalize(antimicrobial),
         class = capitalize(class))

meta_filt <- meta %>% 
  filter(high_microbe_count) %>%
  filter(hap_vap_cap %in% c("HAP", "VAP"))

df <- fread("results/amr_out/amr_matrices/resfinder.amr.decontam.agg.csv") %>%
  filter(run_id %in% meta_filt$run_id) %>%
  pivot_longer(!run_id, names_to = "antimicrobial", values_to = "presence") %>%
  inner_join(meta_filt %>% select(run_id, hap_vap_cap, hap_vap2)) %>%
  mutate(antimicrobial = capitalize(antimicrobial))

meta_filt %>%
  filter(!(run_id %in% df$run_id))

n_total <- n_distinct(df$run_id)

plot_df <- df %>%
  group_by(antimicrobial) %>%
  summarise(n = sum(presence == 1)) %>%
  arrange(desc(n)) %>%
  mutate(prop = n / n_total) %>%
  left_join(drug_meta) %>%
  mutate(class = case_when(antimicrobial == "Streptogramin_b" ~ "Streptogramin",
                           antimicrobial == "Streptogramin_a" ~ "Streptogramin",
                           antimicrobial == "Clindamycin/lincomycin" ~ "Lincosamide",
                           antimicrobial == "Aminopenicillin" ~ "Beta-lactam",
                           antimicrobial == "Amoxicillin_clavulanic acid" ~ "Beta-lactam",
                           TRUE ~ class))

order_df <- plot_df %>%
  head(20)

pal <- c("#BFCDD9FF", "#486078FF","#E9D097FF", "grey50","#AE93BEFF", "pink", "#3FB8AFFF", "#4A7169FF", "#49271BFF")

plot_df %>% 
  filter(antimicrobial %in% order_df$antimicrobial) %>%
  mutate(antimicrobial = factor(antimicrobial, unique(order_df$antimicrobial))) %>%
  mutate(class = factor(class, unique(plot_df$class))) %>%
  ggplot(aes(x = antimicrobial, y = prop, fill = class)) +
  geom_col(color = "black") +
  scale_fill_manual(values = pal, na.value = "white") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Resistance", y = "Prop. patients")

ggsave("results/amr_out/most_common_antibiotics.pdf", dpi = 600, width = 12, height = 3)

# Calculate prop. with >=1 resistance
mat <- fread("results/amr_out/amr_matrices/resfinder.amr.decontam.csv") %>%
  filter(run_id %in% meta_filt$run_id)

row_sums <- rowSums(mat %>% column_to_rownames("run_id")) 

sum(row_sums != 0) / length(row_sums)
tibble(run_id = names(row_sums), n_classes = row_sums) %>%
  left_join(meta_filt) %>%
  ggplot(aes(x = hap_vap2, y = n_classes)) +
  geom_boxplot()

row_sums %>%
  ggplot(aes(x = ))

plot_df %>% View()
