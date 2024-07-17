rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(doParallel)
require(randomcoloR)

meta <- fread("data/metadata/parsed_patient_metadata.filt.csv")

check_df <- fread("results/assembly_out/checkm2_out/quality_report.tsv") %>%
  rename_all(~tolower(.x)) %>%
  dplyr::rename(genome_name = name)

tax_df <- fread("results/assembly_out/tax_out/classify/gtdbtk.bac120.summary.tsv") %>%
  dplyr::rename(genome_name = user_genome)

merged <- check_df %>%
  filter(completeness >= 90, 
         contamination <= 5) %>%
  # Parse taxonomy
  left_join(tax_df) %>%
  mutate(classification = gsub("d__|p__|c__|o__|f__|g__|s__", "", classification)) %>%
  separate(classification, 
           c("domain", "phylum", "class", 
             "order", "family", "genus", 
             "species"), "\\;") %>%
  # Parse run names
  separate(genome_name, c("run", "barcode"), "-", remove = F) %>%
  mutate(run = gsub("Library|INHALE_FRESH_", "", run),
         barcode = gsub("barcode|barcode0", "", barcode)) %>%
  mutate(barcode = gsub("a", "", barcode)) %>%
  mutate(run_id = str_glue("{run}_{barcode}")) %>%
  left_join(meta %>% select(run_id, hospital, hap_vap_cap)) %>%
  select(genome_name, run_id, hospital, hap_vap_cap, 
         family, genus, species, 
         closest_genome_reference, closest_genome_ani,
         checkm2_completeness = completeness, 
         checkm2_contamination = contamination,
         completeness_model_used, genome_size,
         gc_content, coding_density, contig_n50, 
         average_gene_length, total_coding_sequences,
         checkm2_notes = additional_notes, gtdbtk_warnings = warnings) %>%
  # Remove ambiguous taxonomic assignments
  filter(gtdbtk_warnings == "N/A") %>%
  mutate(hospital = ifelse(grepl("Library", genome_name), "Local GP (controls)", hospital))

# Stratify by species and site
plot_df <- merged %>%
  group_by(species, hospital) %>%
  summarise(n = n()) %>%
  arrange(desc(n))

merged %>%
  group_by(hap_vap_cap) %>%
  summarise(n = n())
# Total counts per species
sp_count <- merged %>%
  group_by(species) %>%
  summarise(n = n()) %>%
  arrange(desc(n))

# pal <- distinctColorPalette(n_distinct(plot_df$hospital))
# pal <- c("#A44AE2", "#D4C7CC", "#DF777A", "#84D9BE", "#87A4DB", "#DBCA7B", "#9BE565", "#CE76CA")
pal <- list(UCLH = "darkcyan", RFH = "tan1", Cromwell = "darkorchid4", 
                GSTT = "#2e4057", `Local GP (controls)` = "khaki", GOSH = "cornflowerblue",
                `North Mid` = "#d1495b", `ChelWest` = "cadetblue2")

plot_df %>%
  mutate(species = factor(species, unique(sp_count$species))) %>%
  ggplot(aes(x = species, y = n, fill = hospital)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = pal) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top") +
  labs(x = "Species", y = "No. assemblies", fill = "Site")

# Genus level
plot_genus_df <- merged %>%
  group_by(genus, hospital) %>%
  summarise(n = n()) %>%
  arrange(desc(n))

genus_count <- merged %>%
  group_by(genus) %>%
  summarise(n = n()) %>%
  arrange(desc(n))

# pal <- distinctColorPalette(n_distinct(plot_genus_df$hospital))

plot_genus_df %>%
  mutate(genus = factor(genus, unique(genus_count$genus))) %>%
  ggplot(aes(x = genus, y = n, fill = hospital)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = pal) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top") +
  labs(x = "Species", y = "No. assemblies", fill = "Site")

ggsave("results/assembly_out/assembly_summary_counts.pdf", dpi = 600, width = 7, height = 3)

merged %>%
  group_by(genus) %>%
  summarise(n = n()) %>%
  mutate(genus = factor(genus, unique(genus_count$genus))) %>%
  ggplot(aes(x = genus, y = n, fill = genus)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = distinctColorPalette(n_distinct(genus_count$genus))) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(x = "Species", y = "High-quality completeassemblies", fill = "Site")

# Plot quality
check_df %>% 
  ggplot(aes(x = completeness)) +
  geom_histogram()
# Extract genome list for top pathogens
merged %>%
  filter(species == "Staphylococcus aureus") %>%
  # filter(species == "Klebsiella pneumoniae") %>%
  mutate(genome_path = str_glue("{genome_name}.fna")) %>%
  select(genome_path) %>%
  fwrite("data/metadata/bug_metadata/s_aureus.assemblies.txt",
         eol = "\n",
         col.names = F)

merged %>%
  filter(species == "Escherichia coli") %>%
  # filter(species == "Moraxella catarrhalis") %>%
  # filter(species == "Klebsiella pneumoniae") %>%
  mutate(genome_path = str_glue("{genome_name}.fna")) %>%
  select(genome_path) %>%
  fwrite("data/metadata/bug_metadata/e_coli.assemblies.txt",
         eol = "\n",
         col.names = F)

merged %>%
  summarise(median_n50 = median(contig_n50),
            lower_n50 = quantile(contig_n50, 0.25),
            upper_n50 = quantile(contig_n50, 0.75))
