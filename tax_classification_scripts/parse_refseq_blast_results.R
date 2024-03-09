rm(list = ls())
setwd("C:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)

meta <- fread("data/metadata/parsed_patient_metadata.csv")
# Parse taxonomy metadata
tax_meta <- fread("data/metadata/refseq_genomes_151223/refseq_genomes_151223.metadata.tsv") %>%
  rename_all(~tolower(gsub(" ", "_", .x))) %>%
  filter(assembly_accession != "Assembly Accession") %>%
  select(assembly_accession, 
         ani_best_ani_match_organism, 
         ani_best_ani_match_ani,
         organism_name) %>%
  mutate(ani_best_ani_match_ani = as.numeric(ani_best_ani_match_ani)) %>%
  filter(ani_best_ani_match_ani > 95 | is.na(ani_best_ani_match_ani)) %>%
  mutate(organism_name = ifelse(is.na(ani_best_ani_match_ani), 
                                organism_name,
                                ani_best_ani_match_organism)) %>%
  select(-ani_best_ani_match_ani, -ani_best_ani_match_organism)


tax_meta %>% nrow() == tax_meta %>% distinct(assembly_accession) %>% nrow()

# Assembly accession to chromosome hookup
acc2chrom <- fread("data/metadata/refseq_genomes_151223/ref_to_chrom_mapping.csv",
                   header = F) %>%
  rename(assembly_accession = V1, chrom = V2) %>%
  filter(assembly_accession %in% tax_meta$assembly_accession)

res_dir <- "results/tax_classification_out/blast_out/"

file_list <- list.files(res_dir, full.names = T)

set.seed(66)

morsels <- foreach(file_name = file_list) %do% {
  # file_name = file_list[6]
  # Check if file size is zero
  if (file.size(file_name) != 0) {
    run <- gsub(res_dir, "", file_name)
    run <- gsub("INHALE_FRESH_|barcode0|barcode|.tsv", "", run)
    run <- gsub("a|A", "", run, ignore.case = T)
    run <- gsub("-", "_", run)
  
    temp <- fread(file_name) %>%
      mutate(run_id = run) %>%
      as_tibble()
    
    colnames(temp) <- c("contig_id", "chrom", "pident", 
                       "aln_len", "mismatch", "gapopen", 
                       "qstart", "qend", "sstart", 
                       "send", "evalue", "bitscore",
                       "run_id")
    
    # Retain only hits where reference taxon is resolved
    temp <- temp %>%
      filter(chrom %in% unique(acc2chrom$chrom))
    
    # Get best hit for each contig
    iter_df <- temp %>%
      distinct(run_id, contig_id)
    
    iter_morsels <- foreach(i = seq(nrow(iter_df))) %do% {
      row <- iter_df[i, ]
      best_temp <- temp %>%
        filter(run_id == row$run_id, 
               contig_id == row$contig_id)
      
      max_bitscore <- max(best_temp$bitscore)
      
      # get only top hit
      best_filt <- best_temp %>%
        filter(bitscore == max_bitscore) %>%
        # Get first entry
        head(1)
        
        # Randomly select if more than one hit with same bitscore
        # sample_n(1, replace = F)
      
      return(best_filt)
    }
    
    return(bind_rows(iter_morsels))
  } else {
    return(NULL)
  }
}

dat <- bind_rows(morsels)

parsed <- dat %>%
  left_join(acc2chrom) %>%
  left_join(tax_meta) %>%
  relocate(run_id, .before = 1)

parsed %>% 
  left_join(meta) %>%
  ggplot(aes(x = organism_name, y = run_id, fill = pident)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3)) +
  facet_grid(rows = vars(hap_vap_cap), space = "free", scale = "free")

ggsave("results/test.pdf", width = 30, height = 30)
fwrite(parsed, "results/tax_classification_out/parsed_refseq_blast.best_hits.tsv")


parsed %>%
  inner_join(meta) %>%
  filter(hap_vap_cap == "Water control") %>%
  filter(organism_name == "Escherichia coli") %>% 
  distinct(assembly_accession) %>%View()
