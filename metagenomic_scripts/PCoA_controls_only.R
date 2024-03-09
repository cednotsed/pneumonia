setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(vegan)
require(ape)

meta <- fread("data/metadata/parsed_patient_metadata.csv")

controls <- meta %>%
  filter(hap_vap_cap == "Water control")

RA_filt <- fread("results/metagenomic_out/RA.S.zeroed.csv") %>%
  filter(run_id %in% controls$run_id) %>%
  column_to_rownames("run_id")

as.matrix(RA_filt)[1:5, 1:5]
nrow(RA_filt)

# Bray-Curtis PCoA
bc <- vegdist(as.matrix(RA_filt), method = "bray")

PCOA <- pcoa(bc)

eigenvals <- round(PCOA$values$Relative_eig * 100, 1)

plot_df <- as.data.frame(PCOA$vectors) %>%
  rownames_to_column("run_id") %>%
  left_join(meta) %>%
  mutate(run_date = as.Date(run_date, "%d/%m/%Y"))

# Format dates
brk <- seq.Date(min(plot_df$run_date), max(plot_df$run_date), by = "3 month")

plot_df %>%
  ggplot(aes(Axis.1, Axis.2, color = run_date)) +
  geom_point(size = 5) +
  # scale_color_manual(values = c("indianred", "blue", "goldenrod", "grey50")) +
  theme_classic() + 
  labs(x = str_glue("PCo1 ({eigenvals[1]}%)"),
       y = str_glue("PCo2 ({eigenvals[2]}%)"),
       color = "Run date",
       shape = "Run date",
       title = "Bray-Curtis PCoA") +
  scale_color_gradient2(trans = "date", midpoint = as.numeric(mean(plot_df$run_date)), 
                        breaks = brk,
                        labels = format(brk, "%Y-%b")) 
ggsave("results/metagenomic_out/PCoA.controls_only.png",
       height = 5, width = 5)

pca <- prcomp(RA_filt, retx = T)

as.data.frame(pca$x) %>%
  rownames_to_column("run_id") %>%
  left_join(meta) %>%
  ggplot(aes(x = PC1, y = PC3, color = as.Date(run_date))) +
  geom_point()

