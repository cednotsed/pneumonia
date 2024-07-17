rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(doParallel)

meta <- fread("data/metadata/parsed_patient_metadata.csv")

dir_path <- "results/tax_classification_out/kraken2_out/"
file_list <- list.files(dir_path, full.names = T)
file_list <- file_list[grepl(".tsv", file_list)]
id_list <- gsub(".tsv", "", list.files(dir_path))
head(file_list)
head(id_list)

rank <- "S"

# Find all unique taxa
taxa_df <- fread("databases/k2_pluspfp_20240112/inspect.txt")

# Remove plants
taxa_df %>%
  filter(V4 == "K")

plant_start <- which(taxa_df$V6 == "Viridiplantae")
plant_end <- which(taxa_df$V6 == "Metazoa")

to_remove <- taxa_df[plant_start:(plant_end - 1), ] %>%
  distinct(V6)

taxa_filt <- deframe(taxa_df %>% 
                       filter(V4 == rank) %>%
                       distinct(V6) %>%
                       filter(!(V6 %in% to_remove$V6)))

taxa_list <- c("unclassified", taxa_filt)

# cl <- makeCluster(12)
# registerDoParallel(cl)
# 
# # Time function
start <- Sys.time()

morsels <- foreach(i = seq(length(file_list)), 
                   .packages = c("tidyverse", "data.table")) %do% {
  # i = 13
  # rank <- "G"
  file_path <- file_list[i]
  id <- gsub(dir_path, "", file_path)
  id <- gsub(".tsv", "", id)
  
  # Create taxa_list
  temp_list <- c(rep(0, length(taxa_list) + 1))
  names(temp_list) <- c("run", taxa_list)
  
  # Add id
  temp_list[["run"]] <- id
  
  # Read kraken report
  temp <- fread(file_list[i])
    
  # Account for empty reports
  if (nrow(temp) > 0) {
    temp_filt <- temp %>%
      filter(V6 %in% c("U", rank)) %>%
      filter(V8 %in% taxa_list)
    
    # RECORD READS!!!
    temp_list[temp_filt$V8] <- as.numeric(temp_filt$V2)
  } 
  
  return(as_tibble(t(as.data.frame(temp_list))))
}

# stopCluster(cl)
end <- Sys.time()

print(end - start)

final <- bind_rows(morsels)

# Parse col names
# colnames(final) <- gsub("[[:punct:]]", "", colnames(final))
# colnames(final) <- tolower(gsub(" ", "_", colnames(final)))

# Remove zero columns
final_parsed <- final %>%
  column_to_rownames("run") %>%
  mutate(across(everything(), as.numeric))

to_keep <- colSums(final_parsed) != 0

final_parsed <- final_parsed[, to_keep] %>%
  rownames_to_column("run") %>%
  # Parse run ids
  filter(!grepl("EXTRACTION_TEST", run)) %>%
  mutate(run = gsub("\\.no_human|Library|_raw|INHALE_FRESH_|barcode0|barcode", "", run)) %>%
  mutate(run = gsub("a|A", "", run, ignore.case = T)) %>%
  mutate(run_id = gsub("-", "_", run)) %>%
  select(-run) %>%
  # Remove runs not in meta
  filter(run_id %in% meta$run_id)
  
fwrite(final_parsed, str_glue("results/tax_classification_out/abundance_matrices/abundance_matrix.{rank}.tsv"),
       row.names = F)

dim(final_parsed)
