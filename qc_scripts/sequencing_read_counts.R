setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(Hmisc)

meta <- fread("data/metadata/parsed_patient_metadata.csv")

file_dir <- "data/sequencing_summaries/"
file_list <- list.files(file_dir, full.names = T)

set.seed(66)
file_list <- sample(file_list, 10, replace = F)

read_df <- foreach(file = file_list, .combine = "bind_rows") %do% {
  # file = file_list[1]
  id <- gsub(file_dir, "", file)
  id <- gsub("_sequencing_summary.parsed.csv.gz", "", id)
  
  temp <- fread(file) %>%
    filter(barcode_arrangement != "unclassified") %>%
    mutate(run = id)
  
  return(temp)
}

read_df %>%
  ggplot(aes(x = mean_qscore_template)) +
  geom_density() +
  labs(x = "Mean quality", y = "Read density") +
  theme_classic()
# # PASS FAIL UNCLASSIFIED
# plot_df <- read_df %>%
#   mutate(is_unclassified = barcode_arrangement == "unclassified") %>%
#   mutate(passes_filtering = ifelse(passes_filtering, "Pass", "Fail")) %>%
#   mutate(read_type = ifelse(is_unclassified, 
#                             "Unclassified", 
#                             capitalize(passes_filtering))) %>%
#   group_by(run, read_type) %>%
#   summarise(total_reads = sum(n_reads)) %>%
#   mutate(run = gsub("INHALE_FRESH_", "", run)) %>%
#   arrange(desc(total_reads))
# 
# plot_df %>%
#   mutate(run = factor(run, unique(plot_df$run))) %>%
#   mutate(read_type = factor(read_type, c("Unclassified", "Fail", "Pass"))) %>%
#   ggplot(aes(x = run, y = total_reads, fill = read_type)) +
#   geom_bar(stat = "identity") +
#   scale_fill_manual(values = c("grey40", "salmon", "steelblue"))+
#   theme_bw() +
#   labs(x = "Run", y = "Read count", fill = "Read classification") 
# 
# plot_df %>%
#   pivot_wider(id_cols = run, 
#               names_from = "read_type", 
#               values_from = "total_reads") %>%
#   mutate(perc = Pass / sum(Pass + Fail + Unclassified) * 100) %>%
#   arrange(perc) %>%
#   ggplot(aes(perc)) +
#   geom_histogram(bins = 50) +
#   labs(x = "% reads passed filtering (Q>8)", y = "No. runs")

# By barcode
plot_df2 <- read_df %>% 
  filter(barcode_arrangement != "unclassified") %>%
  mutate(barcode_arrangement = gsub("barcode|a|A", "", barcode_arrangement)) %>%
  mutate(barcode_arrangement = as.numeric(barcode_arrangement)) %>%
  mutate(run = gsub("Library|_raw|INHALE_FRESH_", "", run)) %>% 
  mutate(run_id = str_glue("{run}_{barcode_arrangement}")) %>% 
  filter(run_id %in% meta$run_id) %>%
  filter(passes_filtering) %>%
  group_by(run_id) %>%
  summarise(n_reads = n_distinct(read_id),
            median_read_length = median(sequence_length_template),
            sdev_read_length = sd(sequence_length_template))

plot_df2 %>% 
  fwrite("results/qc_out/sample_read_counts.csv")

plot_df2 %>%
  filter(n_reads >= 1000) %>%
  select(run_id) %>%
  fwrite("results/qc_out/high_read_samples.csv")

## Plots ##
plot_df2 %>%
  ggplot(aes(x = log10(n_reads))) +
  geom_histogram() +
  geom_vline(xintercept = 3,
             color = "red",
             lty = "dashed") +
  labs(x = "Log10(read count)", y = "No. samples")

# Qscore per sample
plot_df2 %>%
  ggplot(aes(x = median_qscore)) +
  geom_histogram() +
  labs(x ="Median Q-score", y = "No. samples")

# Read length per sample
plot_df2 %>%
  ggplot(aes(x = median_length)) +
  geom_histogram() +
  labs(x ="Median read length", y = "No. samples")

# Get summary stats
# Number of samples with >= 1000 reads
plot_df2 %>%
  nrow()

plot_df2 %>%
  filter(n_reads >= 1000) %>%
  nrow()




