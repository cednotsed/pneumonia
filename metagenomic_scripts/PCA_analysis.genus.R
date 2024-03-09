rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(doParallel)
require(randomcoloR)

## FUNCTIONS ##
otu_to_RA <- function(df) {
  mat <- as.matrix(df)
  RA_df <- as.data.frame(mat / rowSums(mat))
  # RA_df <- add_column(RA_df, df$run, .before = 1)
  colnames(RA_df) <- colnames(df)
  
  return(RA_df)
}

otu_to_PA <- function(df, read_threshold) {
  prev_read <- df
  
  prev_read[prev_read <= read_threshold] <- 0
  prev_read[prev_read > read_threshold] <- 1
  colnames(prev_read) <- colnames(df)
  return(prev_read)
}

RA_to_PA <- function(RA_df, PA_threshold) {
  prev_RA <- RA_df
  prev_RA[prev_RA <= PA_threshold] <- 0
  prev_RA[prev_RA > PA_threshold] <- 1
  colnames(prev_RA) <- colnames(RA_df)
  return(prev_RA)
}

meta <- fread("data/metadata/patient_metadata.csv") %>%
  rename_all(~tolower(gsub(" |/", "_", .x))) %>%
  mutate(run = as.character(run)) %>%
  mutate(biofire_organism = ifelse(biofire_valid != "Valid", biofire_organism, "Invalid test")) %>%
  mutate(curetis_organism = ifelse(curetis_organism != "Valid", curetis_organism, "Invalid test")) %>%
  mutate(run_id = str_glue("{run}_{barcode}")) %>%
  mutate(hap_vap_cap = ifelse(grepl("Water", sample_id, ignore.case = T),
                              "Water control",
                              hap_vap_cap)) %>%
  filter(!(sample_id %in% c("Special", "Failed Run - 26")))

table(meta$hap_vap_cap)

dat <- fread("results/metagenomic_out/abundance_matrix.G.tsv") %>%
  select(-`Homo`, -unclassified) %>%
  filter(!grepl("EXTRACTION_TEST", run)) %>%
  as_tibble() %>%
  mutate(run = gsub("INHALE_FRESH_|barcode0|barcode", "", run)) %>%
  mutate(run_id = gsub("-", "_", run)) %>%
  select(-run) %>%
  # Remove runs not in meta
  filter(run_id %in% meta$run_id) %>%
  column_to_rownames("run_id") 

# Filter read counts <= 10, RA <= 0.005
RA_df <- otu_to_RA(dat)
prev_RA <- RA_to_PA(RA_df, 0.005)
prev_read <- otu_to_PA(dat, 10)
prev_df <- as.data.frame(prev_read & prev_RA)

# Set read counts to zero
dat_filt <- dat
dat_filt[!as.matrix(prev_df)] <- 0

# Remove zero taxa
dat_filt <- dat_filt[, colSums(dat_filt) != 0]

# Get relative abundance
RA_filt <- otu_to_RA(dat_filt) %>%
  # Remove IDs with zero reads
  filter(!is.na(`Pseudomonas`)) %>%
  rownames_to_column("run_id") %>%
  # filter(!(run_id %in% c("31_4", "22_8"))) %>%
  column_to_rownames("run_id")

RA_filt[1:5, 1:5]
pca <- prcomp(RA_filt, retx = T)

as.data.frame(pca$x) %>%
  rownames_to_column("run_id") %>%
  left_join(meta) %>%
  ggplot(aes(x = PC1, y = PC3, color = hap_vap_cap)) +
  geom_point()
# geom_text(aes(label = run_id))

long_df <- RA_filt %>%
  rownames_to_column("run_id") %>%
  pivot_longer(!run_id, names_to = "genus", values_to = "rel_a") %>%
  left_join(meta)

summary_stats <- long_df %>%
  group_by(genus) %>%
  summarise(max_abundance = max(rel_a)) %>%
  arrange(desc(max_abundance)) %>%
  head(20)

pal <- distinctColorPalette(n_distinct(summary_stats$genus))
long_df %>%
  filter(genus %in% summary_stats$genus) %>%
  mutate(genus = factor(genus, unique(summary_stats$genus))) %>%
  ggplot(aes(x = run_id, y = rel_a, fill = genus)) +
  geom_bar(stat = "identity", 
           position = "stack",
           color = "black") +
  facet_grid(~hap_vap_cap, scales = "free", space = "free") +
  theme_classic() +
  scale_fill_manual(values = pal) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# merged_df %>% 
#   filter(micro_organism == "Morganella morganii") %>%
#   filter(grepl("Morganella", species))
# 
# RA_df %>%
#   rownames_to_column("run") %>%
#   filter(run == "INHALE_FRESH_1-barcode01") %>%
#   select(`Morganella morganii`)
# filter()
# 
# 
# meta %>%
#   distinct(curetis_valid)
# dat %>%
#   pivot_longer(!c("run", "barcode"), 
#                names_to = "species", 
#                values_to = "read_count") %>%
#   filter(species == "Homo sapiens") %>%
#   arrange(desc(read_count))
# meta %>%
#   select(run, barcode, )
# meta %>%
#   distinct(run)
# colnames(meta)

# # Parse IDs
# RA_parsed <- RA_filt %>%
#   # Remove IDs with zero reads
#   filter(is.na(`Pseudomonas sp. T1-3-2`))
#   # rownames_to_column("run_id")
#   # separate(run, into = c("run", "barcode"), "-") %>%
#   # mutate(barcode = as.numeric(gsub("barcode", "", barcode))) %>%
#   # mutate(barcode = as.character(barcode)) %>%
#   # mutate(barcode = replace_na(barcode, "12a")) %>%
#   # mutate(run = gsub("INHALE_FRESH_", "", run)) %>%
#   # mutate(run_id = str_glue("{run}_{barcode}")) %>%
#   # select(-run, -barcode) %>%
#   # relocate(run_id, .before = 1) %>%
#   # column_to_rownames("run_id")