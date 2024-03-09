setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(doParallel)

dir_path <- "results/metagenomic_out/kraken2_out/"
file_list <- list.files(dir_path, full.names = T)
id_list <- gsub(".tsv", "", list.files(dir_path))
head(file_list)
head(id_list)

rank <- "F"
# Find all unique taxa
taxa_df <- fread("databases/k2_pluspf_20231009/inspect.txt",
                 skip = 7)

taxa_filt <- deframe(taxa_df %>% 
                       filter(V4 == rank) %>%
                       distinct(V6))
taxa_list <- c("unclassified", taxa_filt)

cl <- makeCluster(12)
registerDoParallel(cl)

# Time function
start <- Sys.time()

morsels <- foreach(i = seq(length(file_list)), 
                   .packages = c("tidyverse", "data.table")) %dopar% {
  # i = 13
  rank <- "F"
  file_path <- file_list[i]
  id <- gsub(dir_path, "", file_path)
  id <- gsub(".tsv", "", id)
  
  temp <- fread(file_list[i])
    
  # Create taxa_list
  temp_list <- c(rep(0, length(taxa_list) + 1))
  names(temp_list) <- c("run", taxa_list)
  
  # Add id
  temp_list[["run"]] <- id
  # temp_list[["Pseudomonas aeruginosa"]]
  # Account for empty reports
  if (nrow(temp) > 0) {
    temp_filt <- temp %>%
      filter(V4 %in% c("U", rank)) 
    
    # RECORD READS!!!
    temp_list[temp_filt$V6] <- temp_filt$V2
  } 
  
  return(as_tibble(t(as.data.frame(temp_list))))
}

stopCluster(cl)
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
  rownames_to_column("run")
  

fwrite(final_parsed, str_glue("results/metagenomic_out/abundance_matrix.{rank}.tsv"),
       row.names = F)

dim(final)
# final[1:5, 1:5]
# unique(final["SRR10098296", 1:10])
# fread(paste0(dir_path, "SRR10098296.tsv"))
