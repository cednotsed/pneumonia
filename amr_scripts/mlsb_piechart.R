rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(ggpubr)
require(viridis)
require(scales)
require(ggsci)
require(ggforce)
require(Hmisc)

meta <- fread("data/metadata/parsed_patient_metadata.filt.csv") %>%
  filter(hap_vap_cap %in% c("HAP", "VAP"))

class_meta <- fread("results/amr_out/resfinder_out/resfinder.parsed.decontam.csv") %>%
  distinct(antimicrobial, class) %>%
  mutate(antimicrobial = capitalize(antimicrobial)) %>%
  mutate(class = capitalize(class))

genome_meta <- fread("results/assembly_out/assembly_merged_metadata.csv") %>%
  filter(gtdbtk_warnings == "N/A") %>%
  filter(run_id %in% meta$run_id) %>%
  filter(checkm2_contamination <= 5) %>%
  filter(species != "")

file_dir <- "results/amr_out/resfinder_out.assemblies/all_reports/"
file_list <- list.files(file_dir, full.names = T)

morsels <- foreach(file_name = file_list) %do% {
  id <- gsub(file_dir, "", file_name)
  id <- gsub(".txt", "", id)
  fread(file_name) %>%
    mutate(across(everything(), as.character)) %>%
    mutate(genome_name = id)
}

merged <- bind_rows(morsels) %>%
  select(genome_name, 
         antimicrobial = `# Antimicrobial`, 
         class = Class, 
         wgs_predicted_phenotype = `WGS-predicted phenotype`,
         genetic_background = `Genetic background`) %>%
  mutate(presence = wgs_predicted_phenotype == "Resistant") %>%
  inner_join(genome_meta) %>%
  mutate(antimicrobial = capitalize(antimicrobial),
         class = capitalize(class)) 

plot_df <- merged %>%
  filter(wgs_predicted_phenotype == "Resistant") %>%
  filter(class %in% c("Macrolide", "Streptogramin b", "Lincosamide")) %>%
  separate(genetic_background, c("gene"), sep = "\\(") %>% 
  mutate(gene = ifelse(gene %in% c("erm", "msr"), "msr/erm", gene)) %>%
  group_by(gene) %>%
  summarise(n = n_distinct(run_id)) %>%
  ungroup() %>%
  mutate(total = sum(n))

rpie <- 1
rlabel <-  0.6 * rpie

plot_df %>%
  mutate(end_angle = 2 * pi * cumsum(n) / total,      # ending angle for each pie slice
         start_angle = lag(end_angle, default = 0),   # starting angle for each pie slice
         mid_angle = 0.5 * (start_angle + end_angle)) %>%  # middle of each pie slice, for the text label
  ggplot() +
  geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = rpie,
                   start = start_angle, end = end_angle, fill = gene)) +
  geom_text(aes(x = rlabel * sin(mid_angle), y = rlabel * cos(mid_angle), label = round(n, 4)),
            hjust = 0.5, vjust = 1, size = 20) +
  coord_fixed() +
  scale_x_continuous(limits = c(-1, 1), name = "", breaks = NULL, labels = NULL) +
  scale_y_continuous(limits = c(-1, 1), name = "", breaks = NULL, labels = NULL) +
  theme(panel.background = element_rect(fill = "white"),
        plot.margin=grid::unit(c(0,0,0,0),"cm"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks=element_blank(),
        axis.ticks.length = unit(0, "pt")) +
  scale_fill_manual(values = c("#DABD61FF", "deepskyblue4", "grey")) +
  labs(fill = "MLSb genes")

ggsave("results/amr_out/mlsb_piechart.pdf", width = 5, height = 4, dpi = 600)
