rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(vegan)
require(hillR)
require(ggpubr)

tax_meta <- fread("databases/k2_pluspfp_20240112/inspect.txt")
fungi <- tax_meta[501:801, ]$V6
bacteria <- tax_meta[996:20619, ]$V6

meta <- fread("data/metadata/parsed_patient_metadata.filt.csv") %>%
  filter(hap_vap_cap %in% c("HAP", "VAP"))

# Remove fungi
read_filt <- fread("results/tax_classification_out/abundance_matrices/read_counts.G.zeroed.decontam.2.csv") %>%
  filter(run_id %in% meta$run_id) %>%
  column_to_rownames("run_id") %>%
  dplyr::select(any_of(bacteria))

# Rescale relative abundance
otu_to_RA <- function(df) {
  mat <- as.matrix(df)
  RA_df <- as.data.frame(mat / rowSums(mat))
  RA_df[is.na(RA_df)] <- 0
  colnames(RA_df) <- colnames(df)
  
  return(RA_df)
}

bact_RA <- otu_to_RA(read_filt)
bact_RA <- bact_RA[ , colSums(bact_RA) != 0]

fungal_RA <- fread("results/tax_classification_out/abundance_matrices/RA.G.zeroed.decontam.2.csv") %>%
  filter(run_id %in% meta$run_id) %>%
  pivot_longer(!run_id, names_to = "taxa", values_to = "rel_a") %>%
  filter(taxa %in% fungi) %>%
  group_by(run_id) %>%
  summarise(fungal_rel_a = sum(rel_a))

hill_shannon <- hill_taxa(comm = bact_RA, q = 1)

hshan_df <- tibble(run_id = names(hill_shannon),
                   diversity = hill_shannon) %>%
  left_join(meta %>% select(run_id, sample_type, recent_abx)) %>%
  left_join(fungal_RA) %>%
  mutate(recent_abx = gsub("sample", "sampling", recent_abx)) %>%
  mutate(recent_abx = factor(recent_abx, c(NA, "After sampling", "Before sampling")))

plot_df <- hshan_df %>%
  filter(fungal_rel_a > 0)

# Spearman correlation
spearman <- cor.test(plot_df$diversity, plot_df$fungal_rel_a, 
                     method = "spearman") 
rho <- signif(spearman$estimate, 2)
pval <- signif(spearman$p.value, 2)

# Gamma GLM
gammareg <- glm(diversity ~ sample_type + recent_abx + fungal_rel_a,
                data = plot_df,
                family = Gamma())

gamma_sum <- summary(gammareg, dispersion = MASS::gamma.dispersion(gammareg))
gamma_p <- signif(coefficients(gamma_sum)["fungal_rel_a", "Pr(>|z|)"], 2)

plot_df %>%
  ggplot(aes(x = diversity, y = fungal_rel_a)) +
  geom_smooth(method = "glm", 
              method.args = list(family = Gamma()),
              color = "black", fill = "orange") +
  geom_point(color = "steelblue4", alpha = 0.6) +
  theme_bw() +
  labs(x = "Bacterial Hill-Shannon diversity",
       y = "Fungal relative abundance") +
  annotate("text", x = 7, y = 0.8, label = str_glue("Spearman's Rho={rho}, p={pval}\nGamma GLM, p={gamma_p}"))

ggsave("results/metagenomic_out/fungal_versus_bacterial.pdf", dpi = 600, 
       width = 3, height = 3)

