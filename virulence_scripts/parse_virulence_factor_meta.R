rm(list = ls())
setwd("C:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)

gene_meta <- fread("databases/VF_database_050123/VFDB_setA_nt.headers.fna",
                   sep = "\t",
                   header = F) 
meta <- fread("data/metadata/parsed_patient_metadata.filt.csv")

parsed_meta <- gene_meta %>%
  separate(V1, c("part1", "part2", "part3"), "\\[") %>% 
  separate(part1, c("part11", NA, "part13"), "\\(") %>%
  separate(part13, c("part131", "part132"), "\\)") %>% 
  separate(part2, c("part21", "part22", "part23"), "\\(") %>%
  separate(part22, c("part221", "part222"), "\\) - ") %>% 
  mutate(part11 = gsub(">", "", part11)) %>% 
  rename(seq_id = part11,
         vf_id = part221,
         gene = part131,
         gene_full = part132,
         gene_family = part21,
         mechanism = part222,
         mechanism_id = part23,
         organism_name = part3) %>% 
  mutate(mechanism_id = gsub(")]", "", mechanism_id),
         organism_name = gsub("]", "", organism_name))

res_dir <- "results/virulence_out/mapping_out/"

file_list <- list.files(res_dir, full.names = T)
file_list <- file_list[grepl("cov_stat", file_list)]

morsels <- foreach(file_name = file_list) %do% {
  # file_name = file_list[6]
  id <- gsub(res_dir, "", file_name)
  id <- gsub(".cov_stats.tsv.gz", "", id)
  run <- gsub("INHALE_FRESH_|barcode0|barcode", "", id)
  run <- gsub("a|A", "", run, ignore.case = T)
  run <- gsub("-", "_", run)
  
  temp <- fread(file_name) %>%
    filter(coverage > 0) %>%
    mutate(run_id = run)
  
  return(temp)
}

dat <- bind_rows(morsels)

mat <- dat %>%
  filter(coverage > 50) %>%
  select(seq_id = `#rname`, coverage, run_id) %>%
  pivot_wider(id_cols = "run_id", 
              names_from = "seq_id", 
              values_from = "coverage") %>%
  column_to_rownames("run_id")
  

mat[is.na(mat)] <- 0
mat[mat > 0] <- 1

pca <- prcomp(mat, retx = T)

plot_df <- as.data.frame(pca$x) %>%
  rownames_to_column("run_id") %>%
  right_join(meta)

plot_df %>%
  ggplot(aes(x = PC1, y = PC2, color = hap_vap_cap)) +
  geom_point()
       
plot_df %>%
  filter(hap_vap_cap != "Water control") %>%
  ggplot(aes(y = PC1, x = log10(vent_lenght_hours))) +
  geom_point() +
  geom_smooth()

plot_df2 <- dat %>%
  filter(coverage > 50) %>%
  select(seq_id = `#rname`, coverage, run_id) %>%
  right_join(meta) %>%
  left_join(parsed_meta) %>%
  group_by(run_id, n_reads, hap_vap_cap, gene_family, hosp_los_hours) %>%
  summarise(n = n_distinct(seq_id))
  
plot_df2 %>%  
  ggplot(aes(x = gene_family, y = run_id, fill = n)) +
  geom_tile() +
  facet_grid(rows = vars(hap_vap_cap),
             space = "free", 
             scale = "free")


mat2 <- plot_df2 %>%
  ungroup() %>%
  select(run_id, gene_family, n) %>%
  pivot_wider(id_cols = "run_id", names_from = "gene_family", values_from = "n") %>%
  column_to_rownames("run_id")

mat2[is.na(mat2)] <- 0

pca <- prcomp(mat2, retx = T)

pca_df <- as.data.frame(pca$x) %>%
  rownames_to_column("run_id") %>%
  right_join(meta)

pca_df %>%
  ggplot(aes(x = PC1, y = PC3, color = hap_vap_cap)) +
  geom_point()

water_controls <- deframe(meta %>%
  filter(hap_vap_cap == "Water control") %>%
  select(run_id))

control_genes <- dat %>%
  filter(coverage > 50) %>%
  select(seq_id = `#rname`, coverage, run_id) %>%
  filter(run_id %in% water_controls) %>%
  distinct(seq_id)

meta_filt <- meta %>%
  filter(hap_vap_cap != "Water control") %>%
  mutate(hosp_los_hours = as.numeric(hosp_los_hours)) %>%
  filter(!is.na(hosp_los_hours))

mat_filt <- mat[rownames(mat) %in% meta_filt$run_id, 
                !(colnames(mat) %in% control_genes$seq_id)] %>%
  rownames_to_column("run_id") %>%
  left_join(meta_filt)

morsels <- foreach(i = 2:1003) %do% {
  gene_temp <- mat_filt[, i]
  test <- cor.test(gene_temp, as.numeric(mat_filt$hosp_los_hours),
                   na.rm = T,
                   method = "spearman")
  
  tibble(seq_id = colnames(mat_filt)[i],
         p_val = test$p.value,
         correlation = test$estimate)
}

bind_rows(morsels) %>%
  filter(p_val < 0.05) %>%
  arrange(p_val)

mat_filt %>%
  ggplot(aes(x = VFG011424, y = hosp_los_hours)) +
  geom_point()

mat_filt %>%
  filter(VFG011424 == 1) %>%
  distinct(run_id)
mat_filt2 <- mat[, !(colnames(mat) %in% control_genes$seq_id)]
pca <- prcomp(mat_filt2, retx = T)

pca_df <- as.data.frame(pca$x) %>%
  rownames_to_column("run_id") %>%
  right_join(meta_filt)

pca_df %>%
  ggplot(aes(x = PC1, y = PC3, color = log10(hosp_los_hours))) +
  geom_point()
