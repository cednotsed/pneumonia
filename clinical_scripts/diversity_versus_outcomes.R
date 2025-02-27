rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(vegan)
require(hillR)
require(ggpubr)
require(MKpower)

clin_meta <- fread("data/metadata/parsed_clinical_metadata.csv")
  # mutate(parsed_outcome = ifelse(parsed_outcome %in% c("Not cured", "Died"),
  #                                "Not cured or died",
  #                                parsed_outcome))

meta <- fread("data/metadata/parsed_patient_metadata.filt.csv") %>%
  filter(hap_vap_cap %in% c("HAP", "VAP")) %>%
  filter(high_microbe_count)

RA_filt <- fread("results/tax_classification_out/abundance_matrices/RA.G.zeroed.decontam.2.csv") %>%
  filter(run_id %in% meta$run_id) %>%
  column_to_rownames("run_id")

# Remove patients with no microbes
RA_filt <- RA_filt[rowSums(RA_filt) != 0, ]

hill_shan <- hill_taxa(comm = RA_filt, q = 1)

hshan_df <- tibble(run_id = names(hill_shan),
                   diversity = hill_shan)

# Outcomes
clin_df <-  hshan_df %>%
  inner_join(clin_meta) %>%
  mutate(recent_abx = gsub("sample", "sampling", recent_abx)) %>%
  mutate(recent_abx = factor(recent_abx, c(NA, "After sampling", "Before sampling"))) %>%
  filter(!(hosp_los_hours %in% c("", "Not Known"))) %>%
  mutate(hosp_los_hours = as.numeric(hosp_los_hours)) %>%
  mutate(PELOD2 = as.numeric(PELOD2),
         SOFA = as.numeric(SOFA))

outcome_counts <- clin_df %>%
  filter(parsed_outcome != "Unknown") %>%
  group_by(parsed_outcome) %>%
  summarise(n = n())

clin_df %>%
  filter(parsed_outcome != "Unknown") %>%
  ggplot(aes(x = parsed_outcome, y = diversity, fill = parsed_outcome)) +
  geom_boxplot() +
  geom_pwc() +
  geom_text(aes(label = str_glue("n={n}"), y = 9),
            data = outcome_counts) + 
  theme_classic() +
  scale_fill_manual(values = c("royalblue", "grey80", "red")) +
  labs(x = "Outcome at day 21", y = "Hill-Shannon diversity") +
  theme(legend.position = "none")

# ggsave("results/clinical_out/outcome_v_diversity.pdf", dpi = 600, width = 4, height = 4)

# # Power calculation
# hshan_df %>%
#   filter(parsed_outcome != "Unknown") %>%
#   group_by(parsed_outcome) %>%
#   summarise(sd = sd(diversity))
# 
# rx <- function(n) rnorm(n, mean = 0, sd = 2.1) 
# ry <- function(n) rnorm(n, mean = 2, sd = 2.1) 
# ## two-sample
# sim.ssize.wilcox.test(rx = rx, ry = ry, n.max = 100, n.min = 3, step.size = 1, iter = 1000)

# SOFA/PELOD
plot_df <- clin_df %>% 
  mutate(score = ifelse(age_group == "Child", PELOD2, SOFA)) %>%
  filter(!is.na(score)) %>% 
  mutate(score_type = ifelse(age_group == "Child", "PELOD2", "SOFA"))

age_counts <- plot_df %>%
  group_by(score_type) %>%
  summarise(n = n())

plot_df %>%
  ggplot(aes(x = score, y = diversity, fill = score_type)) +
  geom_smooth(aes(color = score_type), fill = "grey") +
  geom_point(pch = 21, color = "black") +
  geom_text(aes(label = str_glue("n={n}"), color = score_type, y = 8, x = 8),
            position = position_dodge(width = 3),
            data = age_counts) +
  scale_color_manual(values = c("darkseagreen4", "hotpink"),
                     guide = "none") +
  scale_fill_manual(values = c("darkseagreen4", "hotpink")) +
  labs(x = "Score", y = "Hill-Shannon diversity", fill = "Organ dysfunction score") +
  theme_classic() +
  theme(legend.position = "none")

ggsave("results/clinical_out/sofa_pelod_diversity.pdf", dpi = 600, width = 5, height = 3)

# Duration
duration_df <- hshan_df %>%
  inner_join(meta) %>%
  mutate(hosp_los_hours = as.numeric(hosp_los_hours),
         vent_length_hours = as.numeric(vent_length_hours))

duration_df %>%
  filter(!is.na(hosp_los_hours)) %>%
  ggplot(aes(x = log10(hosp_los_hours), y = diversity)) +
  geom_smooth(color = "steelblue") +
  geom_point(pch = 21, fill = "steelblue", color = "black", alpha = 0.5) +
  theme_classic() +
  labs(x = "Log10(length of stay, hours)", y = "Hill-Shannon diversity")

ggsave("results/clinical_out/LOS_diversity.pdf", dpi = 600, width = 5, height = 3)

duration_df %>%
  filter(vent_length_hours != 0) %>% 
  ggplot(aes(x = log10(vent_length_hours), y = diversity)) +
  geom_smooth(color = "indianred") +
  geom_point(pch = 21, fill = "indianred", color = "black", alpha = 0.5) +
  theme_classic() +
  labs(x = "Log10(ventilation duration, hours)", y = "Hill-Shannon diversity")

ggsave("results/clinical_out/ventilation_hours_diversity.pdf", dpi = 600, width = 5, height = 3)

# ggpubr::ggarrange(p1, p2, p3, align = "hv", nrow = 1)
# ggsave("results/clinical_out/duration_versus_diversity.pdf", dpi = 600, 
#        height = 3, width = 8)

duration_df %>%
  filter(vent_lenght_hours != 0) %>% 
  nrow()
  
g1 <- glm(diversity ~ hosp_los_hours,
          data = duration_df,
          family = Gamma(link = "identity"))
summary(g1)
# Antibiotics
hshan_df %>%
  ggplot(aes(x = antibiotics, y = diversity)) +
  geom_boxplot() +
  geom_pwc()
greg <- glm(diversity ~ hosp_los_hours + parsed_outcome,
    family = Gamma,
    data = hsimp_df)

# Save results for machine learning
test <- duration_df %>%
  filter(!is.na(hosp_los_hours)) %>%
  filter(hosp_los_hours != 0) %>%
  select(run_id, hosp_los_hours) %>%
  left_join(RA_filt %>% rownames_to_column("run_id")) %>%
  column_to_rownames("run_id") 

to_remove <- apply(test %>% select(-hosp_los_hours), 2, function(x) {sum(x > 0)}) 
to_remove <- names(to_remove)[to_remove < 10]

test %>%
  select(-any_of(to_remove)) %>%
  fwrite("results/test.csv")

l1 <- lm(hosp_los_hours ~ Streptococcus + Rothia + Escherichia + Pseudomonas + Candida + Staphylococcus + Klebsiella + Paraburkholderia + Chryseobacterium + Enterococcus,
   data = test)

summary(l1)
hist(l1$residuals, breaks = 500)
test %>%
  # filter(Streptococcus != 0) %>%
  ggplot(aes(x = Streptococcus, y = log10(hosp_los_hours))) + 
  geom_point() +
  geom_smooth()
summary(l1)
adults <- hsimp_df %>%
  filter(age_group == "Adult") %>%
  mutate(SOFA_class = case_when(SOFA <= 6 ~ "<= 5",
                                # SOFA >= 7 & SOFA <= 9 ~ "7-9",
                                SOFA > 6 ~ ">5")) %>%
  mutate(SOFA_class = factor(SOFA_class, c("<= 5", ">5")))

adults %>%
  ggplot(aes(x = SOFA, y = diversity)) +
  geom_point()


test <- glm(diversity ~ sample_type + antibiotics + SOFA_class,
   data = adults,
   family = Gamma)
summary(test)

adults %>%
  ggplot(aes(x = SOFA_class, y = diversity)) +
  geom_boxplot()
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
