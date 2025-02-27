rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)

hookup <- fread("data/genomes/fungal_references/fungal_references.csv") %>%
  right_join(fread("data/genomes/fungal_references/chrom_to_acc_hookup.csv"))

df <- fread("results/tax_classification_out/abundance_matrices/read_counts.S.zeroed.csv") %>%
  pivot_longer(!run_id, names_to = "species", values_to = "k2_reads")

file_dir <- "results/fungal_out/stats_out.comp/"
file_list <- list.files(file_dir, "tsv", full.names = T)

morsels <- foreach(file_name = file_list) %do% {
  fread(file_name) %>%
    mutate(file_name = file_name)
}

merged <- bind_rows(morsels) %>%
  mutate(file_name = gsub(str_glue("{file_dir}|.tsv"), "", file_name)) %>%
  separate(file_name, c("run", "barcode"), "-") %>%
  mutate(run = gsub("INHALE_FRESH_", "", run),
         barcode = gsub("barcode|barcode0|a", "", barcode)) %>%
  mutate(run_id = str_glue("{run}_{barcode}")) %>%
  rename(chrom = `#rname`)

plot_df <- merged %>%
  left_join(hookup) %>%
  group_by(run_id, species) %>%
  summarise(mapped_reads = sum(numreads),
            coverage = sum(covbases) / sum(endpos)) %>%
  inner_join(df)
  
fungal_df <- fread("results/fungal_out/fungal_runs.csv")   

cor.test(plot_df$mapped_reads, plot_df$k2_reads)

plot_df %>%
  ggplot(aes(x = mapped_reads, y = k2_reads)) +
    geom_point() 

plot_df %>%
  ungroup() %>%
  distinct(species)
plot_df %>%
  filter(mapped_reads > 10000, 
         k2_reads == 0)

plot_df %>%
  ggplot(aes(x = coverage, y = log10(mapped_reads + 1))) +
  geom_point()

plot_df %>%
  filter(coverage > 0.5) %>% 
  distinct(run_id) %>% View()
taxa <- unique(fungal_df$taxa)
taxa_list <- paste0(taxa, collapse = "|")

df %>%
  filter(run_id %in% fungal_df$run_id) %>%
  filter(k2_reads != 0) %>%
  filter(grepl(taxa_list, species)) %>%
  distinct(species)  %>% View()

fungal_df %>% arrange(desc(rel_a)) %>% View()

df <-fread("results/fungal_out/fungal_runs.csv")
df %>%
  filter(grepl("Candida|Naka|Asper", taxa)) %>%
  distinct(run_id)
