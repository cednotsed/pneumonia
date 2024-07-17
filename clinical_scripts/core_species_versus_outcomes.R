rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(vegan)
require(hillR)
require(ggpubr)

meta <- fread("data/metadata/parsed_clinical_metadata.csv") %>%
  left_join(fread("data/metadata/antibiotic_metadata.csv")) %>%
  filter(parsed_outcome != "Unknown")
# mutate(parsed_outcome = ifelse(parsed_outcome %in% c("Not cured", "Died"), 
#                                "Not cured or died",
#                                parsed_outcome))

RA_filt <- fread("results/tax_classification_out/abundance_matrices/RA.G.zeroed.decontam.2.csv") %>%
  filter(run_id %in% meta$run_id) %>%
  select(run_id, Streptococcus, Rothia, Prevotella, Fusobacterium) %>%
  left_join(meta) %>%
  mutate(recent_abx = gsub("sample", "sampling", recent_abx)) %>%
  mutate(recent_abx = factor(recent_abx, c(NA, "After sampling", "Before sampling"))) %>%
  filter(!(hosp_los_hours %in% c("", "Not Known"))) %>%
  mutate(hosp_los_hours = as.numeric(hosp_los_hours)) %>%
  mutate(PELOD2 = as.numeric(PELOD2),
         SOFA = as.numeric(SOFA))
  

linreg <- lm(Streptococcus ~ sample_type + hosp_los_hours + vent_lenght_hours + hap_vap_cap + recent_abx + PELOD2,
   data = RA_filt %>%
     filter(!is.na(PELOD2)))
summary(linreg)
RA_filt %>%
  ggplot(aes(x = Prevotella, y = log10(PELOD2))) +
  geom_point() +
  geom_smooth()

linreg <- lm(PELOD2 ~ recent_abx + Streptococcus + Rothia + Prevotella + Fusobacterium,
   data = RA_filt)

summary(linreg)
summary(linreg)
left_join(meta) %>%
  mutate(recent_abx = gsub("sample", "sampling", recent_abx)) %>%
  mutate(recent_abx = factor(recent_abx, c(NA, "After sampling", "Before sampling"))) %>%
  filter(!(hosp_los_hours %in% c("", "Not Known"))) %>%
  mutate(hosp_los_hours = as.numeric(hosp_los_hours)) %>%
  mutate(PELOD2 = as.numeric(PELOD2),
         SOFA = as.numeric(SOFA))

adults <- hsimp_df %>%
  filter(age_group == "Adult") %>%
  mutate(SOFA_class = case_when(SOFA <= 5 ~ "<= 5",
                                # SOFA >= 7 & SOFA <= 9 ~ "7-9",
                                SOFA > 5 ~ ">5")) %>%
  mutate(SOFA_class = factor(SOFA_class, c("<= 5", ">5")))

outcome_counts <- adults %>%
  group_by(parsed_outcome) %>%
  summarise(n = n())

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
