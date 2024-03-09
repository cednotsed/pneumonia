setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(vegan)
require(ape)

meta <- fread("data/metadata/parsed_patient_metadata.csv")

RA_filt <- fread("results/metagenomic_out/RA.S.zeroed.decontam.csv") %>%
  left_join(meta %>% select(run_id, hap_vap_cap)) %>%
  filter(hap_vap_cap != "Water control") %>%
  select(-hap_vap_cap) %>%
  column_to_rownames("run_id")

RA_filt <- RA_filt[rowSums(RA_filt) != 0, ]



as.matrix(RA_filt)[1:5, 1:5]
dim(RA_filt)
# Bray-Curtis PCoA
bc <- vegdist(as.matrix(RA_filt), method = "bray")

PCOA <- pcoa(bc)

eigenvals <- round(PCOA$values$Relative_eig * 100, 1)

plot_df <- as.data.frame(PCOA$vectors) %>%
  rownames_to_column("run_id") %>%
  left_join(meta)

plot_df %>%
  ggplot(aes(Axis.1, Axis.2, color = hap_vap_cap, shape = hap_vap_cap)) +
  geom_point() +
  scale_color_manual(values = c("indianred", "blue", "goldenrod", "grey50")) +
  theme_classic() + 
  labs(x = str_glue("PCo1 ({eigenvals[1]}%)"),
       y = str_glue("PCo2 ({eigenvals[2]}%)"),
       color = "Sample Type",
       shape = "Sample Type",
       title = "Bray-Curtis PCoA")

ggsave("results/metagenomic_out/PCoA.decontaminated.png",
       height = 5, width = 5)
