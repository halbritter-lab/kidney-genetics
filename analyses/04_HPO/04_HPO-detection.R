############################################
## load libraries
library(tidyverse)  ## needed for general table operations
library(jsonlite)   ## needed for HGNC requests
library(rvest)      ## needed for scraping
library(httr)       ## needed for scraping
library(kableExtra) ## needed to present scraped data
library("R.utils")  ## gzip downloaded and result files
library(config)     ## needed for config loading
############################################


############################################
## define relative script path
project_topic <- "nephrology"
project_name <- "kidney-genetics"
script_path <- "/analyses/04_HPO/"

## read configs
config_vars_proj <- config::get(file = Sys.getenv("CONFIG_FILE"),
    config = project_topic)

## set working directory
setwd(paste0(config_vars_proj$projectsdir, project_name, script_path))

## set global options
options(scipen = 999)
############################################


############################################
# load global functions
# hgnc functions
source("../functions/hgnc-functions.R", local = TRUE)
source("../functions/hpo-functions.R", local = TRUE)
source("../functions/file-functions.R", local = TRUE)
# helper functions
source("../functions/helper-functions.R", local = TRUE)
############################################


############################################
# compute date only once or somehow in config
current_date <- get_current_date_iso8601()
############################################

############################################
## download all required database sources from HPO
# we load and use the results of previous walks through the ontology tree if not older then 1 month

# HPO obo file download
if (check_file_age("hpo_obo", "../shared/data/downloads/", 1)) {
  hpo_obo_filename <- get_newest_file("hpo_obo", "../shared/data/downloads/")
} else {
  # HPO links to hpo_obo file needs to be set in config
  hpo_obo_url <- config_vars_proj$hpo_obo_url

  hpo_obo_filename <- paste0("../shared/data/downloads/hpo_obo.",
    current_date,
    ".obo")

  download.file(hpo_obo_url, hpo_obo_filename, mode = "wb")
}

hpo <- get_ontology(
    hpo_obo_filename,
    propagate_relationships = "is_a",
    extract_tags = "everything",
    merge_equivalent_terms = TRUE
)
############################################


############################################
## get relevant HPO terms for kidney disease classification

# get all children of term upper urinary tract (HP:0010935) for classification into kidney disease groups
# walk through the ontology tree and add all unique terms descending from
# Abnormality of the upper urinary tract (HP:0010935)
# we load and use the results of previous walks through the ontology tree if not older then 1 month

if (check_file_age("hpo_list_kidney", "../shared/", 1)) {
  hpo_list_kidney <- read_csv(get_newest_file("hpo_list_kidney", "../shared"))
} else {
  all_hpo_children_list_kidney <- hpo_all_children_from_term("HP:0010935", hpo)

  # transform the list into a tibble
  hpo_list_kidney <- all_hpo_children_list_kidney %>%
    unlist() %>%
    tibble(`term` = .) %>%
    unique() %>%
    mutate(query_date = get_current_date_iso8601())

  write_csv(hpo_list_kidney,
    file = paste0("../shared/hpo_list_kidney.",
      get_current_date_iso8601(),
      ".csv"),
    na = "NULL")

  gzip(paste0("../shared/hpo_list_kidney.", get_current_date_iso8601(), ".csv"),
    overwrite = TRUE)
}
############################################


############################################
## download all required database sources from HPO and OMIM
# we load and use the results of previous walks through the ontology tree if not older then 1 month

# disease ontology annotations from HPO
if (check_file_age("phenotype", "../shared/data/downloads/", 1)) {
  phenotype_hpoa_filename <- get_newest_file("phenotype", "../shared/data/downloads/")
} else {

  # disease ontology annotations from HPO

  phenotype_hpoa_filename <- paste0("../shared/data/downloads/phenotype.",
    get_current_date_iso8601(),
    ".hpoa")

  download.file(config_vars_proj$phenotype_hpoa_url, phenotype_hpoa_filename, mode = "wb")

  gzip(phenotype_hpoa_filename,
    overwrite = TRUE)

  phenotype_hpoa_filename <- paste0(phenotype_hpoa_filename,
    ".gz")
}

# OMIM links to genemap2 file needs to be set in config and applied for at
# https://www.omim.org/downloads
if (check_file_age("omim_genemap2", "../shared/data/downloads/", 1)) {
  omim_genemap2_filename <- get_newest_file("omim_genemap2", "../shared/data/downloads/")
} else {

  # OMIM links to genemap2 file needs to be set in config and applied for at
  # https://www.omim.org/downloads
  omim_genemap2_url <- config_vars_proj$omim_genemap2_url

  omim_genemap2_filename <- paste0("../shared/data/downloads/omim_genemap2.",
    get_current_date_iso8601(),
    ".txt")

  download.file(omim_genemap2_url, omim_genemap2_filename, mode = "wb")

  gzip(omim_genemap2_filename,
    overwrite = TRUE)

  omim_genemap2_filename <- paste0(omim_genemap2_filename,
    ".gz")
}
############################################


############################################
## search OMIM and Orphanet for HPO terms and filter
phenotype_hpoa <- read_delim(phenotype_hpoa_filename,
    delim = "\t",
  escape_double = FALSE,
    trim_ws = TRUE,
  skip = 4)

omim_genemap2 <- read_delim(omim_genemap2_filename, "\t",
    escape_double = FALSE,
    col_names = FALSE,
    comment = "#",
    trim_ws = TRUE) %>%
  select(Chromosome = X1,
    Genomic_Position_Start = X2,
    Genomic_Position_End = X3,
    Cyto_Location = X4,
    Computed_Cyto_Location = X5,
    MIM_Number = X6,
    Gene_Symbols = X7,
    Gene_Name = X8,
    approved_symbol = X9,
    Entrez_Gene_ID = X10,
    Ensembl_Gene_ID = X11,
    Comments = X12,
    Phenotypes = X13,
    Mouse_Gene_Symbol_ID = X14) %>%
  select(approved_symbol, Phenotypes) %>%
  separate_rows(Phenotypes, sep = "; ") %>%
  separate(Phenotypes, c("disease_ontology_name", "hpo_mode_of_inheritance_term_name"), "\\), (?!.+\\))") %>%
  separate(disease_ontology_name, c("disease_ontology_name", "Mapping_key"), "\\((?!.+\\()") %>%
  mutate(Mapping_key = str_replace_all(Mapping_key, "\\)", "")) %>%
  separate(disease_ontology_name, c("disease_ontology_name", "MIM_Number"), ", (?=[0-9][0-9][0-9][0-9][0-9][0-9])") %>%
  mutate(Mapping_key = str_replace_all(Mapping_key, " ", "")) %>%
  mutate(MIM_Number = str_replace_all(MIM_Number, " ", "")) %>%
  filter(!is.na(MIM_Number))  %>%
  filter(!is.na(approved_symbol)) %>%
  mutate(disease_ontology_id = paste0("OMIM:", MIM_Number)) %>%
  separate_rows(hpo_mode_of_inheritance_term_name, sep = ", ") %>%
  mutate(hpo_mode_of_inheritance_term_name = str_replace_all(hpo_mode_of_inheritance_term_name, "\\?", "")) %>%
  select(-MIM_Number) %>%
  unique() %>%
  mutate(hpo_mode_of_inheritance_term_name = case_when(hpo_mode_of_inheritance_term_name == "Autosomal dominant" ~ "Autosomal dominant inheritance",
    hpo_mode_of_inheritance_term_name == "Autosomal recessive" ~ "Autosomal recessive inheritance",
    hpo_mode_of_inheritance_term_name == "Digenic dominant" ~ "Digenic inheritance",
    hpo_mode_of_inheritance_term_name == "Digenic recessive" ~ "Digenic inheritance",
    hpo_mode_of_inheritance_term_name == "Isolated cases" ~ "Sporadic",
    hpo_mode_of_inheritance_term_name == "Mitochondrial" ~ "Mitochondrial inheritance",
    hpo_mode_of_inheritance_term_name == "Multifactorial" ~ "Multifactorial inheritance",
    hpo_mode_of_inheritance_term_name == "Pseudoautosomal dominant" ~ "X-linked dominant inheritance",
    hpo_mode_of_inheritance_term_name == "Pseudoautosomal recessive" ~ "X-linked recessive inheritance",
    hpo_mode_of_inheritance_term_name == "Somatic mosaicism" ~ "Somatic mosaicism",
    hpo_mode_of_inheritance_term_name == "Somatic mutation" ~ "Somatic mutation",
    hpo_mode_of_inheritance_term_name == "X-linked" ~ "X-linked inheritance",
    hpo_mode_of_inheritance_term_name == "X-linked dominant" ~ "X-linked dominant inheritance",
    hpo_mode_of_inheritance_term_name == "X-linked recessive" ~ "X-linked recessive inheritance",
    hpo_mode_of_inheritance_term_name == "Y-linked" ~ "Y-linked inheritance"))

phenotype_hpoa_filter <- phenotype_hpoa %>%
   filter(hpo_id %in% hpo_list_kidney$term) %>%
   select(database_id, hpo_id) %>%
   unique() %>%
  group_by(database_id) %>%
  summarise(abnormality_of_the_kidney_hpo_terms = paste(hpo_id, collapse = " | "),
    .groups = "keep") %>%
  ungroup()

hpo_gene_list <- phenotype_hpoa_filter %>%
  left_join(omim_genemap2, by = c("database_id" = "disease_ontology_id")) %>%
  filter(!is.na(approved_symbol)) %>%
  mutate(database_and_hpo_id = paste0(database_id, " (", abnormality_of_the_kidney_hpo_terms, ")")) %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(approved_symbol)) %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  select(-disease_ontology_name, -Mapping_key, -hpo_mode_of_inheritance_term_name) %>%
  group_by(approved_symbol) %>%
  summarise(approved_symbol = paste(unique(approved_symbol), collapse = " | "),
    hgnc_id = paste(unique(hgnc_id), collapse = " | "),
    abnormality_of_the_kidney_hpo_terms = paste(abnormality_of_the_kidney_hpo_terms, collapse = "; "),
    database_id = paste(database_id, collapse = "; "),
    database_and_hpo_id = paste(database_and_hpo_id, collapse = "; "),
    source_count = n(),
    .groups = "keep") %>%
  ungroup() %>%
  mutate(gene_name_reported = approved_symbol) %>%
  mutate(source_evidence = (source_count > 0)) %>%
  select(approved_symbol,
    hgnc_id,
    gene_name_reported,
    source = database_and_hpo_id,
    source_count,
    source_evidence)

# normalize source_count to 0/1 as percentiles
hpo_gene_list_normalize <- hpo_gene_list %>%
  normalize_percentile("source_count")
############################################


############################################
## save results
write_csv(hpo_gene_list_normalize,
  file = paste0("results/04_HPO_genes.",
    current_date,
    ".csv"),
  na = "NULL")

gzip(paste0("results/04_HPO_genes.", current_date, ".csv"),
  overwrite = TRUE)
############################################
