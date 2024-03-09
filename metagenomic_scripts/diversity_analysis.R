setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(vegan)
require(hillR)
require(ggpubr)

meta <- fread("data/metadata/parsed_patient_metadata.filt.csv")

RA_filt <- fread("results/metagenomic_out/RA.S.zeroed.csv") %>%
  filter(run_id %in% meta$run_id) %>%
  column_to_rownames("run_id")

# Remove zero taxa
RA_filt <- RA_filt[rowSums(RA_filt) != 0, ]

hill_shannon <- hill_taxa(comm = RA_filt, q = 1)

hshan_df <- tibble(run_id = names(hill_shannon),
                diversity = hill_shannon) %>%
  left_join(meta) %>%
  mutate(recent_abx = gsub("sample", "sampling", recent_abx)) %>%
  mutate(recent_abx = factor(recent_abx, c(NA, "After sampling", "Before sampling")))

count_df <- hshan_df %>%
  group_by(hap_vap_cap, recent_abx) %>%
  summarise(n = n())

hshan_df %>%
  ggplot(aes(x = hap_vap_cap, y = diversity, fill = recent_abx)) +
  geom_boxplot() +
  geom_pwc() +
  geom_text(aes(x = hap_vap_cap, y = 0, label = str_glue("n={n}")),
            data = count_df,
            position = position_dodge(width = 0.8)) +
  labs(x = "Patient group", 
       y = "Hill-Shannon Diversity (Richness)",
       fill = "Antibiotics given") +
  theme_bw() +
  scale_fill_discrete(na.value = "grey")

shannon_lr <- glm((as.numeric(factor(hap_vap_cap)) - 1) ~ recent_abx + diversity,
   data = hshan_df,
   family = "binomial")

summary(shannon_lr)

hill_simpson <- hill_taxa(comm = RA_filt, q = 2)

hsimp_df <- tibble(run_id = names(hill_simpson),
                diversity = hill_simpson) %>%
  left_join(meta)
  # filter(hap_vap_cap != "CAP")

count_df <- hsimp_df %>%
  group_by(hap_vap_cap, recent_abx) %>%
  summarise(n = n())

hsimp_df %>% 
  # mutate(hap_vap_cap = ifelse(hap_vap_cap == "Water control", hap_vap_cap, "Pneumonia")) %>%
  ggplot(aes(x = hap_vap_cap, y = log10(diversity))) +
  geom_violin() +
  geom_pwc()
  geom_text(aes(x = hap_vap_cap, y = 1, label = str_glue("n={n}")),
            data = count_df,
            position = position_dodge(width = 0.8)) +
  labs(x = "Patient group", 
       y = "Hill-Simpson Diversity (Evenness)",
       fill = "Antibiotics given")

hsimp_df %>% 
  filter(is.na(hap_vap_cap))
shannon_lr <- glm((as.numeric(factor(hap_vap_cap)) - 1) ~ recent_abx + diversity,
                  data = hshan_df,
                  family = "binomial")

count_df <- hsimp_df %>%
  group_by(hap_vap_cap, recent_abx, sample_type) %>%
  summarise(n_samples = n()) %>%
  ungroup() %>%
  complete(hap_vap_cap, recent_abx, sample_type) %>%
  mutate(n_samples = replace_na(n_samples, 0))

# hsimp_df %>%
#   ggplot(aes(x = hap_vap_cap, y = diversity, fill = recent_abx)) +
#   geom_boxplot() +
#   facet_grid(.~sample_type) +
#   geom_text(aes(x = hap_vap_cap, 
#                 y = Inf, 
#                 label = str_glue("n={n_samples}"),
#                 color = recent_abx),
#             position = position_dodge(width = 0.9),
#             vjust = 2,
#             size = 2,
#             data = count_df) +
#   geom_pwc()

hsimp_df %>%
  filter(hap_vap_cap %in% c("Water control", "VAP")) %>%
  # mutate(hap_vap_cap = ifelse(hap_vap_cap == "Water control", hap_vap_cap, "Pneumonia")) %>%
  group_by(hap_vap_cap) %>%
  summarise(mean_div = mean(diversity),
            sd = sd(diversity),
            lower_q = quantile(diversity, 0.25),
            upper_q = quantile(diversity, 0.75),
            min = min(diversity),
            max = max(diversity))

test <- hsimp_df %>%
  filter(hap_vap_cap == "Water control")

shapiro.test(test$diversity)
