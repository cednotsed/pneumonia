rm(list = ls())
setwd("C:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)

meta <- fread("data/metadata/parsed_patient_metadata.csv")

meta %>%
  filter(hap_vap_cap == "VAP") %>%
  distinct(sample_id, .keep_all = T) %>%
  ggplot(aes(x = vent_lenght_hours / 24)) +
  geom_density(bins = 50, fill = "grey") +
  geom_vline(aes(xintercept = 0, color = "0"), lty = "dashed") +
  geom_vline(aes(xintercept = 2, color = "2"), lty = "dashed") +
  geom_vline(aes(xintercept = 4, color = "4"), lty = "dashed") +
  geom_vline(aes(xintercept = 6, color = "6"), lty = "dashed") +
  geom_vline(aes(xintercept = 8, color = "8"), lty = "dashed") +
  theme_classic() +
  scale_color_manual(values = c("cadetblue1", "deepskyblue", "dodgerblue2",
                                "dodgerblue4", "blue4")) +
  labs(x = "Duration of ventilation (days)", 
       y = "Patient density",
       color = "Timepoint (days)") +
  theme(legend.position = "right")
  
  
plot_df <- meta %>%
  filter(hap_vap_cap == "VAP") %>%
  distinct(sample_id, .keep_all = T) %>%
  mutate(vent_days = vent_lenght_hours / 24)
  
duration <- plot_df$vent_days
breaks <- seq(2, 30, by=0.1) 
duration.cut <- cut(duration, breaks, right = F)
duration.freq <- table(duration.cut)
  
duration.freq
cumfreq0 <- c(0, cumsum(duration.freq)) 
tibble(x = breaks, y = cumfreq0 / nrow(plot_df)) %>%
  ggplot(aes(x, y)) +
  geom_line() +
  geom_vline(aes(xintercept = 0, color = "0"), lty = "dashed") +
  geom_vline(aes(xintercept = 2, color = "2"), lty = "dashed") +
  geom_vline(aes(xintercept = 4, color = "4"), lty = "dashed") +
  geom_vline(aes(xintercept = 6, color = "6"), lty = "dashed") +
  geom_vline(aes(xintercept = 8, color = "8"), lty = "dashed") +
  theme_classic() +
  scale_color_manual(values = c("cadetblue1", "deepskyblue", "dodgerblue2",
                                "dodgerblue4", "blue4")) +
  labs(x = "Duration of ventilation (days)", 
       y = "Cum. prop. of patients",
       color = "Timepoint (days)") +
  theme(legend.position = "right")

plot_df %>%
  summarise(sum(vent_days <= 8) / n())

quantile(plot_df$vent_days, 0.95)
