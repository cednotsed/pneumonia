rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(segmented)
require(Biostrings)

fna <- readDNAStringSet("data/genomes/all_irep_refs.no_plasmids.fna")
chrom_names <- str_split(names(fna), "\\ ", simplify = T)[, 1]
len_df <- tibble(chromosome = chrom_names, ref_length = width(fna))

cov_df <- fread("results/irep_out/mapping_out/INHALE_FRESH_11-barcode01.cov_stats.txt.gz") %>%
  arrange(V1, V2) %>%
  dplyr::rename(chromosome = V1, pos = V2, coverage = V3)

len_filt <- len_df %>%
  filter(chromosome %in% cov_df$chromosome)

all_pos <- foreach(chrom = unique(len_filt$chromosome), 
                   .combine = "bind_rows") %do% {
  temp <- len_filt %>%
    filter(chromosome == chrom)
  
  tibble(pos = seq(temp$ref_length)) %>%
    mutate(chromosome = chrom)
}
nrow(all_pos)
len_filt
cov_df %>% 
  filter(coverage == 0)

cov_df %>%
  group_by(Chromosome) %>%
  summarise(n = n()) %>%
  arrange(desc(n))

window_size <- 10000
window_step <- 10

# Calculate average coverage for each window
average_coverage <- cov_df %>%
  group_by(chromosome, window = ceiling(row_number() / window_step)) %>%
  summarise(average_coverage = mean(coverage))

# Drop any incomplete windows
average_coverage <- average_coverage[complete.cases(average_coverage), ]
average_coverage

# Filter out windows with low coverage (adjust threshold as needed)
threshold <- 8
filtered_coverage <- average_coverage %>%
  filter(average_coverage > threshold)  # Set a suitable threshold

# Log2-transform coverage values
filtered_coverage$log2_coverage <- log2(filtered_coverage$average_coverage)

# Fit piecewise linear function
seg_model <- segmented(lm(log2_coverage ~ Window, data = filtered_coverage),
                       seg.Z = ~Window, psi = list(Window = c(5000)),
                       control = seg.control(display = FALSE))

# Predict breakpoints
breakpoints <- segmented.lm.breakpoints(seg_model)
origin <- breakpoints[1, "Estimate"]
terminus <- breakpoints[2, "Estimate"]


seg_model <- segmented(lm(log2_coverage ~ Window, 
                          data = filtered_coverage %>%
                            filter(Chromosome == "NZ_CP007265.1")),
                       seg.Z = ~Window, psi = list(Window = c(5000)),
                       control = seg.control(display = FALSE))

filtered_coverage %>%
  filter(chromosome == "NZ_CP082270.1") %>%
  ggplot(aes(x = window, y = log2_coverage)) +
  geom_line() 

filtered_coverage %>%
  ggplot(aes(x = window, y = log2_coverage)) +
  geom_line() 

cov_df %>%
  group_by(chromosome) %>%
  summarise(n = n()) %>%
  arrange(desc(n))
