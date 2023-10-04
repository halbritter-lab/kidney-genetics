############################################
## load libraries
library(tidyverse)  ## needed for general table operations
library(biomaRt)    ## needed to get gene coordinates
library(STRINGdb)   ## needed to compute StringDB identifiers
library("R.utils")  ## gzip downloaded and result files
library(readxl)     ## needed to read xlsx file
library(config)     ## needed for config loading
library(janitor)    ## needed for cleaning column names
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
source("../functions/gtex-functions.R", local = TRUE)
source("../functions/file-functions.R", local = TRUE)
############################################


############################################
## download all required database sources from HGNC, OMIM, gnomAD, clinvar and genCC
# we load and use the results of previous walks through the ontology tree if not older then 1 month

current_date <- strftime(as.POSIXlt(Sys.time(),
    "UTC", "%Y-%m-%dT%H:%M:%S"),
  "%Y-%m-%d")

# HGNC file download
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

  non_alt_loci_set_filename <- paste0(non_alt_loci_set_filename,
    ".gz")
}

# OMIM file download
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

  omim_genemap2_filename <- paste0(omim_genemap2_filename,
    ".gz")
}

# gnomAD lof metrics download
if (check_file_age("gnomad.v2.1.1.lof_metrics.by_gene", "../shared/data/downloads/", 1)) {
  gnomad_v211_lof_met_gene_filename <- get_newest_file("gnomad.v2.1.1.lof_metrics.by_gene", "../shared/data/downloads/")
} else {
  # gnomAD lof metrics links to genemap2 file needs to be set in config
  gnomad_v211_lof_met_gene_url <- config_vars_proj$gnomad_v211_lof_met_gene_url

  gnomad_v211_lof_met_gene_filename <- paste0("../shared/data/downloads/gnomad.v2.1.1.lof_metrics.by_gene.",
    current_date,
    ".txt.gz")

  download.file(gnomad_v211_lof_met_gene_url, gnomad_v211_lof_met_gene_filename, mode = "wb")
}

# clinvar VCF file download
if (check_file_age("clinvar", "../shared/data/downloads/", 1)) {
  clinvar_vcf_hg19_filename <- get_newest_file("clinvar", "../shared/data/downloads/")
} else {
  # clinvar VCF file links to genemap2 file needs to be set in config
  clinvar_vcf_hg19_url <- config_vars_proj$clinvar_vcf_url

  clinvar_vcf_hg19_filename <- paste0("../shared/data/downloads/clinvar.",
    current_date,
    ".vcf.gz")

  download.file(clinvar_vcf_hg19_url, clinvar_vcf_hg19_filename, mode = "wb")
}

# genCC file download
if (check_file_age("gencc_submissions", "../shared/data/downloads/", 1)) {
  gencc_submissions_filename <- get_newest_file("gencc_submissions", "../shared/data/downloads/")
} else {
  # genCC file links to genemap2 file needs to be set in config
  gencc_submissions_url <- config_vars_proj$gencc_submissions_url

  gencc_submissions_filename <- paste0("../shared/data/downloads/gencc_submissions.",
    current_date,
    ".xlsx")

  download.file(gencc_submissions_url, gencc_submissions_filename, mode = "wb")
}

# ClinGen file download
if (check_file_age("clingen_gene_disease_summary", "../shared/data/downloads/", 1)) {
  clingen_gene_disease_summary_filename <- get_newest_file("clingen_gene_disease_summary", "../shared/data/downloads/")
} else {
  # ClinGen file links to genemap2 file needs to be set in config
  clingen_gene_disease_summary_url <- config_vars_proj$clingen_gene_disease_summary_url

  clingen_gene_disease_summary_filename <- paste0("../shared/data/downloads/clingen_gene_disease_summary.",
    current_date,
    ".csv")

  download.file(clingen_gene_disease_summary_url, clingen_gene_disease_summary_filename, mode = "wb")

  gzip(clingen_gene_disease_summary_filename,
    overwrite = TRUE)

  clingen_gene_disease_summary_filename <- paste0(clingen_gene_disease_summary_filename,
    ".gz")
}
############################################


############################################
## load the downloaded HGNC file
# TODO: specify column specifications to suppress warnings
non_alt_loci_set <- read_delim(paste0(non_alt_loci_set_filename),
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
## load gnomAD lof metrics file
gnomad_v211_lof_met_gene <- read_delim(gnomad_v211_lof_met_gene_filename,
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE)
############################################


############################################
# load or process clinvar VCF file into table
# and then safe/load the result as a table for future use
# TODO: maybe safe the clinvar table in the (shared) or analysis data folder instead
if (check_file_age("clinvar_table", "results/", 1)) {
  clinvar_table_filename <- get_newest_file("clinvar_table", "results/")

  clinvar_table <- read_csv(clinvar_table_filename)
} else {
  ## load the downloaded clinvar VCF file
  clinvar_vcf <- read_delim(clinvar_vcf_hg19_filename,
      delim = "\t", escape_double = FALSE,
      trim_ws = TRUE,
      comment = "##")

  # reformat the vcf into a table
  clinvar_table <- clinvar_vcf %>% 
    separate_rows(INFO, sep = ";") %>%
    separate(INFO, into = c("key", "value"), sep = "=") %>%
    spread(key, value) %>%
    separate_rows(GENEINFO, sep = "\\|") %>%
    separate_rows(CLNSIG, sep = "/")

  # write the table to a csv file
  write_csv(clinvar_table,
    file = paste0("results/clinvar_table.", current_date, ".csv"))

  gzip(paste0("results/clinvar_table.", current_date, ".csv"),
    overwrite = TRUE)
}
############################################


############################################
# load and process genCC file
# Define the order of classification_title
classification_order <- c(
  "Definitive",
  "Strong",
  "Moderate",
  "Limited",
  "Supportive",
  "Disputed Evidence",
  "Refuted Evidence",
  "No Known Disease Relationship"
)

# load the Excel file
gencc_submissions_table <-  read_excel(gencc_submissions_filename) %>%
  mutate(classification_title = factor(classification_title, levels = classification_order))
############################################


############################################
# load and process ClinGen CSV file
# clean the column names
# remove comments
clingen_gene_disease_summary_table <- read_csv(clingen_gene_disease_summary_filename, 
    skip = 4) %>%
  janitor::clean_names() %>%
  filter(!str_detect(gene_symbol, "\\+"))
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
# add ensembl_gene_id_version for hg19 and hg38 using the ensembl functions
# add genocode id using the GTEx API and the get_multiple_gencode_ids function
# TODO: fix warning "! Ensembl will soon enforce the use of https. Ensure the 'host' argument includes https://""
non_alt_loci_set_coordinates <- non_alt_loci_set_string %>%
  mutate(ensembl_gene_id_version_hg19 =
    gene_id_version_from_ensembl(ensembl_gene_id, reference = "hg19")$ensembl_gene_id_version) %>%
  mutate(ensembl_gene_id_version_hg38 =
    gene_id_version_from_ensembl(ensembl_gene_id, reference = "hg38")$ensembl_gene_id_version) %>%
  mutate(gencode_id = get_multiple_gencode_ids(symbol)$gencode_id) %>%
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
# annotate gnomAD pLI and missense Z-scores
# currently using: use gnomAD download table
# TODO: future adaption use https://gnomad.broadinstitute.org/api/

# split the mane_select column into mane_enst and mane_nm to find the corresponding transcript for the join
non_alt_loci_set_coordinates_reformat <- non_alt_loci_set_coordinates %>%
  separate(mane_select, into = c("mane_enst", "mane_nm"), sep = "\\|") %>%
  separate(mane_enst, into = c("mane_enst", "mane_enst_version"), sep = "\\.") %>%
  dplyr::select(symbol, mane_enst) %>%
  filter(!is.na(mane_enst))

# first join with the above helper table
gnomad_v211_lof_met_gene_mane_enst <- non_alt_loci_set_coordinates_reformat %>%
  left_join(gnomad_v211_lof_met_gene, by = c("mane_enst" = "transcript")) %>%
  dplyr::select(symbol, mane_enst, pLI, pNull, pRec, oe_lof, oe_lof_lower, oe_lof_upper, oe_mis, oe_syn, lof_z, mis_z, syn_z, constraint_flag) %>%
  mutate_at(vars(-symbol, -mane_enst, -constraint_flag), ~ if_else(is.na(.), NA_real_, round(as.numeric(.), 6)))

# then join with non_alt_loci_set_coordinates_omim
non_alt_loci_set_coordinates_gnomad <- non_alt_loci_set_coordinates %>%
  left_join(gnomad_v211_lof_met_gene_mane_enst, by = c("symbol"))
############################################


############################################
# annotate with OMIM P numbers
# use omim tables

# group the omim table by gene symbol
omim_genemap2_grouped <- omim_genemap2 %>%
  group_by(approved_symbol) %>%
  summarise(
    omim_summary = paste(paste0(disease_ontology_name,
      "[",
      disease_ontology_id,
      "]-",
      hpo_mode_of_inheritance_term_name), collapse = " | "),
    omim_count = n(),
    .groups = 'drop') %>%
  ungroup()

# joining with non_alt_loci_set_coordinates_gnomad
non_alt_loci_set_coordinates_gnomad_omim <- non_alt_loci_set_coordinates_gnomad %>%
  left_join(omim_genemap2_grouped, by = c("symbol" = "approved_symbol"))
############################################


############################################
# annotate with GenCC classification
# first summarize the table by gene symbol
gencc_submissions_table_grouped <- gencc_submissions_table %>%
  arrange(gene_symbol, classification_title) %>%
  group_by(gene_symbol, gene_curie, disease_curie, disease_title, moi_title) %>%
  summarise(
    classification_title = paste0(classification_title, collapse = ", "),
    gencc_classification_count = n(),
    .groups = 'drop') %>%
  ungroup() %>%
  group_by(gene_symbol, gene_curie) %>%
  summarise(
    gencc_summary = paste(paste0(classification_title,
      " (",
      disease_title,
      "[",
      disease_curie,
      "]-",
      moi_title,
      ")"), collapse = " | "),
    gencc_count = paste(gencc_classification_count,
      collapse = ", "),
    .groups = 'drop') %>%
  ungroup() %>%
  dplyr::select(-gene_symbol)

# join with non_alt_loci_set_coordinates_gnomad_omim
non_alt_loci_set_coordinates_gnomad_omim_gencc <- non_alt_loci_set_coordinates_gnomad_omim %>%
  left_join(gencc_submissions_table_grouped, by = c("hgnc_id" = "gene_curie"))
############################################


############################################
# annotate ClinGen Gene-Disease Validity
# first summarize the table by gene symbol
clingen_gene_disease_summary_table_grouped <- clingen_gene_disease_summary_table %>%
  arrange(gene_symbol, classification) %>%
  group_by(gene_symbol, gene_id_hgnc, disease_id_mondo, disease_label, moi) %>%
  summarise(
    classification = paste0(classification, collapse = ", "),
    clingen_classification_count = n(),
    .groups = 'drop') %>%
  ungroup() %>%
  group_by(gene_symbol, gene_id_hgnc) %>%
  summarise(
    clingen_summary = paste(paste0(classification,
      " (",
      disease_label,
      "[",
      disease_id_mondo,
      "]-",
      moi,
      ")"), collapse = " | "),
    clingen_count = paste(clingen_classification_count,
      collapse = ", "),
    .groups = 'drop') %>%
  ungroup() %>%
  dplyr::select(-gene_symbol)

# join with non_alt_loci_set_coordinates_gnomad_omim_gencc
non_alt_loci_set_coordinates_gnomad_omim_gencc_clingen <- non_alt_loci_set_coordinates_gnomad_omim_gencc %>%
  left_join(clingen_gene_disease_summary_table_grouped, by = c("hgnc_id" = "gene_id_hgnc"))

############################################


############################################
# annotate clinvar variant counts
# use the clinvar VCF file

# summarize the clinvar table first
# TODO: implement computing HGNC and symbol using hgnc functions
# TODO: join on something safe then like the HGNC ID
clinvar_table_counts <- clinvar_table %>%
  group_by(GENEINFO, CLNSIG) %>%
  summarise(count = n(), .groups = 'drop') %>%
  filter(CLNSIG %in% c("Pathogenic", "Likely_pathogenic", "Uncertain_significance")) %>%
  pivot_wider(names_from = CLNSIG, values_from = count, values_fill = 0) %>%
  mutate(Pathogenic = paste0("P:", Pathogenic),
         Likely_pathogenic = paste0("LP:", Likely_pathogenic),
         Uncertain_significance = paste0("VUS:", Uncertain_significance)) %>%
  unite("clinvar", Pathogenic:Uncertain_significance, sep = "; ", remove = FALSE, na.rm = TRUE) %>%
  separate(GENEINFO, into = c("symbol", "Gene_ID"), sep = ":") %>%
  group_by(symbol) %>%
  summarise(Gene_ID = paste(Gene_ID, collapse = "|"),
    clinvar = paste(clinvar, collapse = "|"),
    .groups = 'drop') %>%
  dplyr::select(symbol, clinvar)

# join with non_alt_loci_set_coordinates_gnomad
non_alt_loci_set_coordinates_gnomad_omim_gencc_clingen_clinvar <- non_alt_loci_set_coordinates_gnomad_omim_gencc_clingen %>%
  left_join(clinvar_table_counts, by = c("symbol" = "symbol"))
############################################


############################################
## export table as csv with date of creation
creation_date <- strftime(as.POSIXlt(Sys.time(),
    "UTC", "%Y-%m-%dT%H:%M:%S"),
  "%Y-%m-%d")

write_csv(non_alt_loci_set_coordinates_gnomad_omim_gencc_clingen_clinvar,
  file = paste0("results/non_alt_loci_set_coordinates.", creation_date, ".csv"))

gzip(paste0("results/non_alt_loci_set_coordinates.", creation_date, ".csv"),
  overwrite = TRUE)
############################################