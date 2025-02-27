rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)
require(seqinr)
require(ape)

checkm2_df <- fread("results/assembly_out/assembly_merged_metadata.csv")
meta <- fread("data/metadata/parsed_patient_metadata.filt.csv")

morsels <- foreach(bug_prefix = c("s_aureus", "e_coli", "h_influenzae")) %do% {
  aln <- readDNAStringSet(str_glue("data/alignments/{bug_prefix}.core_gene_alignment_filtered.snps.trimmed.aln"))
  aln_filt <- aln[grepl("INHALE", names(aln))]
  aln_filt <- as.DNAbin(aln_filt)
  
  # subst_rate <- case_when(bug_prefix == "e_coli" ~ 2e-6,
  #                         bug_prefix == "h_influenzae" ~ 2e-6,
  #                         bug_prefix == "s_aureus" ~ 2e-6)
    
  subst_rate <- 2e-6
  aln_length <- case_when(bug_prefix == "e_coli" ~ 2409139,
                          bug_prefix == "h_influenzae" ~ 260312,
                          bug_prefix == "s_aureus" ~ 1581645)
    
  pair_df <- as.data.frame(as.matrix(dist.dna(aln_filt, model = "N"))) %>%
    rownames_to_column("seq1") %>%
    pivot_longer(!seq1, names_to = "seq2", values_to = "dist") %>%
    filter(seq1 != seq2) %>%
    left_join(checkm2_df %>% select(seq1 = genome_name, run_id1 = run_id)) %>%
    left_join(checkm2_df %>% select(seq2 = genome_name, run_id2 = run_id)) %>%
    left_join(meta %>% select(run_id1 = run_id,  date1 = sample_date, hosp1 = hospital, subtype1 = hap_vap2)) %>%
    left_join(meta %>% select(run_id2 = run_id,  date2 = sample_date, hosp2 = hospital, subtype2 = hap_vap2)) %>%
    mutate(date1 = as.Date(date1, "%d/%m/%Y"),
           date2 = as.Date(date2, "%d/%m/%Y")) %>%
    mutate(time_diff = abs(difftime(date1, date2, units = "weeks"))) %>%
    mutate(subst_rate = subst_rate) %>%
    mutate(aln_length = aln_length) %>%
    mutate(expected_e6 = as.numeric(aln_length * time_diff * subst_rate / 52),
           expected_e5 = as.numeric(aln_length * time_diff * subst_rate / 52 * 10)) %>% 
    mutate(same_hospital = hosp1 == hosp2) %>%
    mutate(same_subtype = subtype1 == subtype2) %>%
    mutate(bug = bug_prefix) 
  
  return(pair_df)
}

plot_df <- bind_rows(morsels)
  
plot_df %>% 
  ggplot(aes(x = time_diff, y = log10(dist))) +
  geom_point() +
  geom_point(aes(x = time_diff, y = log10(expected_e6)), color = "red") +
  geom_line(aes(x = time_diff, y = log10(expected_e6)), color = "red") +
  geom_point(aes(x = time_diff, y = log10(expected_e5)), color = "blue") +
  geom_line(aes(x = time_diff, y = log10(expected_e5)), color = "blue") +
  labs(x = "Time diff. (weeks)", y = "log10(pairwise SNPs)") +
  theme_bw() +
  facet_grid(rows = vars(bug))
    

ggsave("results/transmission_out/expected_snps.pdf", dpi = 600, width = 6, height = 4)

# linreg <- lm(log10(dist) ~ same_hospital + same_subtype + time_diff, data = plot_df)
# posreg <- glm(dist ~ same_hospital + same_subtype + time_diff, data = plot_df,
#               family = "poisson")
# summary(posreg)




