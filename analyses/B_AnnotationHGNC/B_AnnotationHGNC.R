############################################
## load libraries
library(tidyverse)  ## needed for general table operations
library(biomaRt)    ## needed to get gene coordinates
library(STRINGdb)  ## needed to compute StringDB identifiers
library("R.utils")  ## gzip downloaded and result files
library(config)     ## needed for config loading
############################################


############################################
## define relative script path
project_topic <- "nephrology"
project_name <- "kidney-genetics"
script_path <- "/analyses/B_AnnotationHGNC/"

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
source("../functions/ensembl-functions.R", local = TRUE)
source("../functions/gnomad-functions.R", local = TRUE)
source("../functions/file-functions.R", local = TRUE)
############################################


############################################
## download all required database sources from HGNC and OMIM
# we load and use the results of previous walks through the ontology tree if not older then 1 month

current_date <- strftime(as.POSIXlt(Sys.time(),
    "UTC", "%Y-%m-%dT%H:%M:%S"),
  "%Y-%m-%d")

if (check_file_age("non_alt_loci_set", "../shared/data/downloads/", 1)) {
  non_alt_loci_set_filename <- get_newest_file("non_alt_loci_set", "../shared/data/downloads/")
} else {
  # HGNC links to non_alt_loci_set file needs to be set in config
  non_alt_loci_set_url <- config_vars_proj$non_alt_loci_set_url

  non_alt_loci_set_filename <- paste0("../shared/data/downloads/non_alt_loci_set.",
    current_date,
    ".txt")

  download.file(non_alt_loci_set_url, non_alt_loci_set_filename, mode = "wb")

  gzip(non_alt_loci_set_filename,
    overwrite = TRUE)
}

if (check_file_age("omim_genemap2", "../shared/data/downloads/", 1)) {
  omim_genemap2_filename <- get_newest_file("omim_genemap2", "../shared/data/downloads/")
} else {
  # OMIM links to genemap2 file needs to be set in config and applied for at
  # https://www.omim.org/downloads
  omim_genemap2_url <- config_vars_proj$omim_genemap2_url

  omim_genemap2_filename <- paste0("../shared/data/downloads/omim_genemap2.",
    current_date,
    ".txt")

  download.file(omim_genemap2_url, omim_genemap2_filename, mode = "wb")

  gzip(omim_genemap2_filename,
    overwrite = TRUE)
}
############################################


############################################
## load the downloaded HGNC file
# TODO: specify column specifications to suppress warnings
non_alt_loci_set <- read_delim(paste0(non_alt_loci_set_filename, ".gz"),
    "\t",
    col_names = TRUE) %>%
  mutate(update_date = current_date)
############################################


############################################
## load OMIM file and reformat them
omim_genemap2 <- read_delim(omim_genemap2_filename, "\t",
    escape_double = FALSE,
    col_names = FALSE,
    comment = "#",
    trim_ws = TRUE) %>%
  dplyr::select(Chromosome = X1,
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
  dplyr::select(approved_symbol, Phenotypes) %>%
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
  dplyr::select(-MIM_Number) %>%
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
############################################


############################################
## load STRINGdb database
string_db <- STRINGdb$new(version = "11.5",
  species = 9606,
  score_threshold = 200,
  input_directory = "data/")
############################################


############################################
## map gene symbols to StringDB identifiers
non_alt_loci_set_table <- non_alt_loci_set %>%
  dplyr::select(symbol) %>%
  unique()

non_alt_loci_set_df <- non_alt_loci_set_table %>%
    as.data.frame()

non_alt_loci_set_mapped <- string_db$map(non_alt_loci_set_df, "symbol")
non_alt_loci_set_mapped_tibble <- as_tibble(non_alt_loci_set_mapped) %>%
  filter(!is.na(STRING_id)) %>%
  group_by(symbol) %>%
  summarise(STRING_id = str_c(STRING_id, collapse = ";")) %>%
  ungroup %>%
  unique()

## join with String identifiers
non_alt_loci_set_string <- non_alt_loci_set %>%
  left_join(non_alt_loci_set_mapped_tibble, by = "symbol")
############################################


############################################
## add gene coordinates from ensembl
# TODO: fix warning "! Ensembl will soon enforce the use of https. Ensure the 'host' argument includes https://""
non_alt_loci_set_coordinates <- non_alt_loci_set_string %>%
  mutate(hg19_coordinates_from_ensembl =
    gene_coordinates_from_ensembl(ensembl_gene_id)) %>%
  mutate(hg19_coordinates_from_symbol =
    gene_coordinates_from_symbol(symbol)) %>%
  mutate(hg38_coordinates_from_ensembl =
    gene_coordinates_from_ensembl(ensembl_gene_id, reference = "hg38")) %>%
  mutate(hg38_coordinates_from_symbol =
    gene_coordinates_from_symbol(symbol, reference = "hg38")) %>%
  mutate(bed_hg19 =
    case_when(
      !is.na(hg19_coordinates_from_ensembl$bed_format) ~
        hg19_coordinates_from_ensembl$bed_format,
      is.na(hg19_coordinates_from_ensembl$bed_format) ~
        hg19_coordinates_from_symbol$bed_format,
    )
  ) %>%
  mutate(bed_hg38 =
    case_when(
      !is.na(hg38_coordinates_from_ensembl$bed_format) ~
        hg38_coordinates_from_ensembl$bed_format,
      is.na(hg38_coordinates_from_ensembl$bed_format) ~
        hg38_coordinates_from_symbol$bed_format,
    )
  ) %>%
  dplyr::select(-hg19_coordinates_from_ensembl,
    -hg19_coordinates_from_symbol,
    -hg38_coordinates_from_ensembl,
    -hg38_coordinates_from_symbol)
############################################



############################################
# TODO: annotate with OMIM P numbers
# use omim tables

############################################



############################################
# TODO: annotate gnomAD pLI and missense Z-scores
# use gnomAD download table from https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_transcript.txt.bgz
# TODO: future adaption use https://gnomad.broadinstitute.org/api/
ensembl_gene_id_gnomad <- non_alt_loci_set_coordinates %>%
    filter(!is.na(ensembl_gene_id)) %>%
    dplyr::select(ensembl_gene_id) %>%
    head(30) %>%
    mutate(gene_data = get_multiple_gene_data_from_gnomad(ensembl_gene_id))
############################################



############################################
# TODO: annotate ClinVar variant counts
# use the clinvar API
# use https://gnomad.broadinstitute.org/api/

############################################



############################################
# TODO: annotate with GeneCC classification

############################################


############################################
## export table as csv with date of creation
creation_date <- strftime(as.POSIXlt(Sys.time(),
    "UTC", "%Y-%m-%dT%H:%M:%S"),
  "%Y-%m-%d")

write_csv(non_alt_loci_set_coordinates,
  file = paste0("results/non_alt_loci_set_coordinates.", creation_date, ".csv"))

gzip(paste0("results/non_alt_loci_set_coordinates.", creation_date, ".csv"),
  overwrite = TRUE)
############################################