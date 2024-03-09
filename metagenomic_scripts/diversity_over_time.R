rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(vegan)
require(hillR)
require(ggpubr)

meta <- fread("data/metadata/parsed_patient_metadata.csv")

RA_filt <- fread("results/metagenomic_out/RA.S.zeroed.csv") %>%
  column_to_rownames("run_id")

# Remove patients with no microbes
RA_filt <- RA_filt[rowSums(RA_filt) != 0, ]

hill_shannon <- hill_taxa(comm = RA_filt, q = 1)

hshan_df <- tibble(run_id = names(hill_shannon),
                   diversity = hill_shannon) %>%
  left_join(meta) %>%
  mutate(recent_abx = gsub("sample", "sampling", recent_abx)) %>%
  mutate(recent_abx = factor(recent_abx, c(NA, "After sampling", "Before sampling"))) %>%
  filter(!(hosp_los_hours %in% c("", "Not Known"))) %>%
  mutate(hosp_los_hours = as.numeric(hosp_los_hours)) 

n_samples <- n_distinct(hshan_df$run_id)

g1 <- glm(diversity ~ sample_type + post_pcr_qubit + hosp_los_hours,
          data = hshan_df,
          family = Gamma(link = "identity"))

g1_sum <- summary(g1)

slope1 <- formatC(g1$coefficients[["hosp_los_hours"]], format = "e", digits = 1)
pval1 <- formatC(g1_sum$coefficients["hosp_los_hours", "Pr(>|t|)"], format = "e", digits = 1)

p1 <- hshan_df %>%
  ggplot() +
  geom_point(aes(x = hosp_los_hours, y = diversity, color = post_pcr_qubit)) +
  geom_smooth(aes(x = hosp_los_hours, y = diversity),
              method = "glm",
              method.args = list(family = Gamma(link = 'identity'))) +
  labs(x = "Hospital LOS (hrs)", y = "Hill-Shannon diversity",
       color = "post-PCR DNA conc.",
       title = str_glue("Effect of Hospital LOS={slope1}, p={pval1}\nn={n_samples}")) +
  theme_bw()

# Ventilator
plot2_df <- hshan_df %>%
  filter(hap_vap_cap == "VAP")

n_vap <- n_distinct(plot2_df$run_id)

g2 <- glm(diversity ~ sample_type + post_pcr_qubit + vent_lenght_hours,
          data = plot2_df,
          family = Gamma(link = "identity"))

g2_sum <- summary(g2)

slope2 <- formatC(g2$coefficients[["vent_lenght_hours"]], format = "e", digits = 1)
pval2 <- formatC(g2_sum$coefficients["vent_lenght_hours", "Pr(>|t|)"], format = "e", digits = 1)

g2_link <- Gamma(link = "identity")

p2 <- plot2_df %>%
  ggplot() +
  geom_point(aes(x = vent_lenght_hours, y = diversity, color = post_pcr_qubit)) +
  geom_smooth(aes(x = vent_lenght_hours, y = diversity),
              method = "glm",
              method.args = list(family = Gamma(link = 'identity'))) +
  labs(x = "Time on ventilator (hrs)", y = "Hill-Shannon diversity",
       color = "post-PCR DNA conc.",
       title = str_glue("Effect of vent. duration={slope2}, p={pval2}\nn={n_vap}")) +
  theme_bw()

ggpubr::ggarrange(p1, p2, common.legend = T)
# hill_simpson <- hill_taxa(comm = RA_filt, q = 2)
# 
# hsimp_df <- tibble(run_id = names(hill_simpson),
#                    diversity = hill_simpson) %>%
#   left_join(meta) %>%
#   mutate(recent_abx = gsub("sample", "sampling", recent_abx)) %>%
#   mutate(recent_abx = factor(recent_abx, c(NA, "After sampling", "Before sampling"))) %>%
#   filter(!(hosp_los_hours %in% c("", "Not Known"))) %>%
#   mutate(hosp_los_hours = as.numeric(hosp_los_hours))
# 
# hsimp_df %>%
#   ggplot(aes(x = duration, y = diversity, color = post_pcr_qubit)) +
#   geom_point() +
#   geom_smooth() +
#   labs(x = "Time on ventilator (hrs)", y = "Hill-Simpson diversity")
# 
# hsimp_df %>%
# lr2 <- lm(diversity ~ sample_type + log(post_pcr_qubit) + log(hosp_los_hours),
#          data = hsimp_df)
# 
# g1 <- glm(diversity ~ sample_type + post_pcr_qubit + hosp_los_hours,
#                    data = hsimp_df,
#                    family = Gamma(link = "log"))
# 
# g2 <- glm(diversity ~ sample_type + post_pcr_qubit + hosp_los_hours,
#                    data = hsimp_df,
#                    family = Gamma(link = "identity"))
# 
# g3 <- glm(diversity ~ sample_type + post_pcr_qubit + hosp_los_hours,
#                    data = hsimp_df,
#                    family = Gamma(link = "inverse"))
# 
# anova(g1, g2, g3, test="LRT")
# summary(aod)
# 
# ahist(gamma_model$residuals)
# hist(hsimp_df$hosp_los_hours)
# hsimp_df %>%
#   filter()
# x <- rnorm(1)