rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(vegan)
require(hillR)
require(ggpubr)

meta <- fread("data/metadata/parsed_clinical_metadata.csv") %>%
  filter(parsed_outcome != "Unknown") %>%
  mutate(parsed_outcome = ifelse(parsed_outcome %in% c("Not cured", "Died"),
                                 "Not cured or died",
                                 parsed_outcome))

micro_meta <- fread("data/metadata/parsed_microbiology_results.csv") %>%
  filter(run_id %in% meta$run_id) %>%
  filter(!(bugs %in% c("Invalid test"))) %>%
  mutate(bugs = ifelse(grepl("Enterobacter", bugs), "Enterobacter spp.", bugs)) %>%
  mutate(bugs = ifelse(grepl("Acinetobacter", bugs), "Acinetobacter spp.", bugs)) %>%
  mutate(bugs = ifelse(grepl("Proteus", bugs), "Proteus spp.", bugs)) %>%
  filter(bugs != "Coliform")

bug_count <- micro_meta %>%
  group_by(run_id) %>%
  summarise(n_bugs = n_distinct(bugs))

mat <- meta %>%
  left_join(bug_count) %>% 
  mutate(poly = ifelse(n_bugs > 1, "Poly", "Mono")) %>%
  select(parsed_outcome, poly) %>%
  group_by(parsed_outcome, poly) %>%
  summarise(n = n()) %>%
  pivot_wider(id_cols = parsed_outcome, names_from = poly, values_from = n)  %>%
  column_to_rownames("parsed_outcome")

fisher.test(mat)

meta
  ggplot(aes(x = parsed_outcome, y = n_bugs)) +
  geom_boxplot()
RA_filt <- fread("results/tax_classification_out/abundance_matrices/RA.G.zeroed.decontam.csv") %>%
  filter(run_id %in% meta$run_id) %>%
  column_to_rownames("run_id")

# Remove patients with no microbes
RA_filt <- RA_filt[rowSums(RA_filt) != 0, ]

hill_simpson <- hill_taxa(comm = RA_filt, q = 2)

hsimp_df <- tibble(run_id = names(hill_simpson),
                   diversity = hill_simpson) %>%
  left_join(meta) %>%
  mutate(recent_abx = gsub("sample", "sampling", recent_abx)) %>%
  mutate(recent_abx = factor(recent_abx, c(NA, "After sampling", "Before sampling"))) %>%
  filter(!(hosp_los_hours %in% c("", "Not Known"))) %>%
  mutate(hosp_los_hours = as.numeric(hosp_los_hours)) %>%
  mutate(PELOD2 = as.numeric(PELOD2),
         SOFA = as.numeric(SOFA))

hsimp_df %>%
  ggplot(aes(x = log10(hosp_los_hours), y = diversity)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "glm", method.args = list(family = Gamma()))

hsimp_df %>%
  filter(vent_lenght_hours != 0) %>%
  ggplot(aes(x = vent_lenght_hours, y = diversity)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "glm", method.args = list(family = Gamma()))

adults <- hsimp_df %>%
  filter(age_group == "Adult") %>%
  mutate(SOFA_class = case_when(SOFA <= 5 ~ "<= 5",
                                # SOFA >= 7 & SOFA <= 9 ~ "7-9",
                                SOFA > 5 ~ ">5")) %>%
  mutate(SOFA_class = factor(SOFA_class, c("<= 5", ">5")))

outcome_counts <- adults %>%
  group_by(parsed_outcome) %>%
  summarise(n = n())

adults %>%
  ggplot(aes(x = parsed_outcome, y = diversity)) +
  geom_boxplot()
# Diversity versus SOFA
adult_mod <- glm(diversity ~ recent_abx + vent_lenght_hours + day_21_outcome + SOFA,
                 data = adults,
                 family = Gamma(link = "log"))

adult_summ <- summary(adult_mod)

slope <- signif(coef(adult_summ)["SOFA", "Estimate"], 3)
t_val <- signif(coef(adult_summ)["SOFA", "t value"], 3)
p_val <- signif(coef(adult_summ)["SOFA", "Pr(>|t|)"], 2)

adults %>% 
  ggplot(aes(x = as.integer(SOFA), y = diversity, color = as.integer(SOFA))) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "SOFA score", y = "Hill-Simpson diversity",
       title = str_glue("Gamma GLM: Effect={slope}, t={t_val}, p={p_val}")) +
  geom_smooth(method = "glm", method.args = list(family = "Gamma")) 

# Diversity versus SOFA class
adults %>%
  ggplot(aes(x = SOFA_class, y = diversity, fill = SOFA_class)) +
  geom_boxplot() +
  geom_pwc() +
  scale_fill_manual(values = c("steelblue4", "indianred")) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "SOFA score", y = "Hill-Simpson diversity")

# Diversity versus outcome
adults %>%
  ggplot(aes(x = parsed_outcome, y = diversity, fill = parsed_outcome)) +
  geom_boxplot() +
  geom_point(position = position_dodge(width = 0.8)) +
  geom_pwc() +
  geom_text(aes(x = parsed_outcome, y = 0, label = str_glue("n={n}")),
            data = outcome_counts) +
  scale_fill_manual(values = c("olivedrab", "grey40", "khaki3")) +
  labs(x = "21-day outcome", y = "Hill-Simpson diversity") +
  theme_classic() +
  theme(legend.position = "none")


# PELOD
children <- hsimp_df %>%
  filter(age_group == "Child") %>%
  mutate(PELOD2_class = case_when(PELOD2 <= 6 ~ "<= 6",
                                  PELOD2 > 6 ~ ">6")) %>%
  mutate(PELOD2_class = factor(PELOD2_class, c("<= 6", ">6")))

# Diversity versus SOFA
children_mod <- glm(diversity ~ PELOD2,
                    data = children,
                    family = Gamma(link = "identity"))

children_summ <- summary(children_mod)

slope <- signif(coef(children_summ)["PELOD2", "Estimate"], 3)
t_val <- signif(coef(children_summ)["PELOD2", "t value"], 3)
p_val <- signif(coef(children_summ)["PELOD2", "Pr(>|t|)"], 2)

children %>% 
  ggplot(aes(x = as.integer(PELOD2), y = diversity, color = as.integer(PELOD2))) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "PELOD-2 score", y = "Hill-Simpson diversity",
       title = str_glue("Gamma GLM: Effect={slope}, t={t_val}, p={p_val}")) +
  geom_smooth(method = "glm", method.args = list(family = "Gamma")) 

# Diversity versus PELOD2 class
children %>%
  filter(!is.na(PELOD2_class)) %>%
  ggplot(aes(x = PELOD2_class, y = diversity, fill = PELOD2_class)) +
  geom_boxplot() +
  geom_pwc() +
  scale_fill_manual(values = c("steelblue4", "indianred")) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "PELOD-2 score", y = "Hill-Simpson diversity")

outcome_counts2 <- children %>%
  group_by(parsed_outcome) %>%
  summarise(n = n())

children %>%
  ggplot(aes(x = parsed_outcome, y = diversity, fill = parsed_outcome)) +
  geom_boxplot() +
  geom_point(position = position_dodge(width = 0.8)) +
  geom_pwc() +
  geom_text(aes(x = parsed_outcome, y = 0, label = str_glue("n={n}")),
            data = outcome_counts2) +
  scale_fill_manual(values = c("olivedrab", "grey40", "khaki3")) +
  labs(x = "21-day outcome", y = "Hill-Simpson diversity") +
  theme_classic() +
  theme(legend.position = "none")

cor.test(adults$SOFA, adults$hosp_los_hours, method = "spearman")
