rm(list = ls())
setwd("c:/git_repos/pneumonia/")
require(tidyverse)
require(data.table)
require(foreach)
require(doParallel)

meta <- fread("data/metadata/hunt_et_al_v0.2/ena_metadata.tsv")
tax_meta <- fread("data/metadata/hunt_et_al_v0.2/species_calls.tsv") %>%
  rename(sample_accession = Sample)

qual_meta <- fread("data/metadata/hunt_et_al_v0.2/checkm2.tsv")

# meta %>%
#   head() %>% View()
species <- "Klebsiella pneumoniae"
path <- "k_pneumoniae.uk.txt"

species <- "Moraxella catarrhalis"
path <- "m_catarrhalis.uk.txt"

meta %>%
  # filter(grepl("Human|Homo sapiens", host, ignore.case = T)) %>%
  # filter(grepl("UK|United Kingdom", country)) %>%
  inner_join(tax_meta) %>%
  filter(Species == species) %>%
  select(sample_accession) %>%
  fwrite(str_glue("data/metadata/bug_metadata/{path}"),
         eol = "\n",
         col.name = F)

# species <- "Escherichia coli"
# path <- "e_coli.uk.txt"
# 
parsed <- meta %>%
  filter(grepl("Human|Homo sapiens", host, ignore.case = T)) %>%
  filter(grepl("UK|United Kingdom", country)) %>%
  filter(!(collection_date %in% c("", "Not available"))) %>%
  separate(collection_date, c("year", "month", "day"), "-", remove = F) %>%
  mutate(year = as.numeric(year),
         month = as.numeric(month),
         day = as.numeric(day)) %>%
  mutate(year_parsed = ifelse(day > 1500, day, year),
         day_parsed = ifelse(day > 1500, year, day))

# parsed %>%
#   inner_join(tax_meta) %>%
#   filter(grepl("Human|Homo sapiens", host, ignore.case = T)) %>%
#   filter(grepl("UK|United Kingdom", country)) %>%
#   filter(Species == species) %>%
#   filter(year_parsed >= 2016, year_parsed <= 2019) %>%
#   select(sample_accession) %>%
#   fwrite(str_glue("data/metadata/bug_metadata/{path}"),
#          eol = "\n",
#          col.name = F)
# parsed %>%
#   distinct(collection_date, year_parsed, year, month, day) %>% View()
  # filter(!is.na(year) & !is.na(month) & !is.na(day)) %>%
  # mutate(parsed_date = ifelse(day > 1500, 
  #                             as.character(as.Date(collection_date, "%d-%m-%Y")),
  #                             as.character(as.Date(collection_date, "%Y-%m-%d")))) %>% 
  # mutate(parsed_date = as.Date(parsed_date)) %>%
  # filter(!is.na(parsed_date))

# # Extract sequences by time
# meta %>%
#   filter(grepl("Human|Homo sapiens", host, ignore.case = T)) %>%
#   filter(grepl("UK|United Kingdom", country)) %>%
#   inner_join(tax_meta) %>%
#   filter(collection_date < "2019-01-01" & collection_date > "2016-01-01")
#   # select(sample_accession)
#   fwrite("data/metadata/bug_metadata/s_aureus.uk.txt")
#   
# meta %>%
#   filter(sample_accession == "SAMD00083522")

