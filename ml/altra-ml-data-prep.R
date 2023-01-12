# R code to pre-process data for ML applications in ALTRA project
# Author: Brenton Graham
# Last Edit: 12/21/22

# Import packages
require(tidyverse)
require(readxl)
require(phyloseq)
require(microbiome)
require(stringr)
require(magrittr)

# Operation for 'not in'
`%not in%` <- Negate(`%in%`)

# Get abbreviated taxa names
abbrev_taxa <- function(physeq) {
  short_names <- physeq %>%
    tax_table() %>% as.data.frame() %>%
    rownames_to_column("full_name") %>%
    mutate(lowest_cat = case_when(
      !is.na(Species) ~ Species,
      !is.na(Genus) ~ Genus,
      !is.na(Family) ~ Family,
      !is.na(Order) ~ Order,
      !is.na(Class) ~ Class,
      !is.na(Phylum) ~ Phylum,
      TRUE ~ Domain)) %>%
    #select(lowest_cat) %>% pull()
    mutate(
      lowest_cat = str_replace_all(lowest_cat, " ", "."),
      phylum_lowest = case_when(
        Phylum == lowest_cat ~ Phylum,
        TRUE ~ paste(Phylum, lowest_cat, sep = "_")),
      phylum_lowest = tolower(phylum_lowest)) %>% 
    select(phylum_lowest) %>% pull()
  return(short_names)
}

# Import data

data <- paste(getwd(), "/../data/ALTRA_clinical&16S.merged.19-Dec-2022.analysis_samples.csv", sep="") %>%
  read_delim(delim = ",") %>%
  mutate(rowname = ifelse(sample_type_16S == "Stool", paste(rowname, "F", sep = ""), rowname)) %>%
  column_to_rownames("rowname")

# OTU table
otu_table <- data %>% select(contains("Bacteria")) %>%
  t() %>% otu_table(taxa_are_rows = T)

# Sample metadata
metadata <- data %>% select(-contains("Bacteria")) %>%
  select(sample_id, ccp3_group, ccp3, everything())

# Taxonomy matrix
otu_names <- otu_table %>% as.data.frame() %>% rownames()
tax_levels <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
otu_tax_split <- sapply(otu_names, FUN=function(x) strsplit(x, "/"))
tax_table <- plyr::ldply(otu_tax_split, rbind)[-1] %>%
  set_rownames(otu_names) %>%
  set_colnames(tax_levels) %>%
  as.matrix()

# Phyloseq object
physeq <- phyloseq(otu_table(otu_table), tax_table(tax_table), sample_data(metadata))

# Output data frame for Python implementation
rel_ab.output <- meta(physeq) %>% select(ccp3_group, ccp3, sample_type_16S, age, gender) %>%
  merge(physeq %>% microbiome::transform("compositional") %>% t() %>% otu_table() %>% as.data.frame(), by = 'row.names') %>%
  filter(sample_type_16S == "Stool") %>%
  filter(ccp3_group != "PosRA") %>%
  select(-sample_type_16S, -ccp3_group, -Row.names) %>%
  mutate(gender = case_when(gender == "Female" ~ 0, gender == "Male" ~ 1))
output_colnames <- abbrev_taxa(physeq) 

# DeepMicro files
## Features (Relative abundance)
rel_ab.output %>%
  select(contains("Bacteria")) %>%
  write.table(paste(getwd(), "/data/altra.stool.deepMicro_features.csv", sep=""), sep = ",", row.names = F, col.names = F)
## Targets
rel_ab.output %>%
  select(ccp3) %>%
  write.table(paste(getwd(), "/data/altra.stool.deepMicro_targets.csv", sep=""), sep = ",", row.names = F, col.names = F)

# Everything else
# Files with targets
rel_ab.output %>%
  set_colnames(c("ccp3", "age", "sex", output_colnames)) %>%
  write.table(paste(getwd(), "/data/altra.stool.rel_ab.csv", sep = ","), row.names = F, col.names = T)
