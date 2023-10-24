############################################
## load libraries
library(readr)      ## needed to read files
library(tidyverse)  ## needed for general table operations
library(rvest)      ## needed for scraping
library(jsonlite)   ## needed for HGNC requests
library("R.utils")  ## gzip downloaded and result files
library(config)     ## needed for config loading
############################################


############################################
## define relative script path
project_topic <- "nephrology"
project_name <- "kidney-genetics"
script_path <- "/analyses/01_PanelApp/"

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
source("../functions/panelapp-functions.R", local = TRUE)
############################################


############################################
## get a list of all PanelApp panels and filter

## PanelApp UK
# store all pages in a list first
baseurl <- "https://panelapp.genomicsengland.co.uk/api/v1/panels/?format=json&page="

# find all pages through API call
panelapp_uk_last_page <- find_last_page(baseurl)

panelapp_uk_pages <- list()
for (i in 1:panelapp_uk_last_page){
  panelapp_uk_page <- fromJSON(paste0(baseurl, i))
  message("Retrieving page ", i)
  panelapp_uk_pages[[i + 1]] <- panelapp_uk_page$results
}

# combine all into one
panelapp_uk <- rbind_pages(panelapp_uk_pages) %>%
  tibble() %>%
  select(id, name, version, version_created) %>%
  mutate(api_call = paste0(
    "https://panelapp.genomicsengland.co.uk/api/v1/panels/",
    id,
    "/?format=json")) %>%
  mutate(panel_source = "panelapp_uk")

## PanelApp Australia
# store all pages in a list first
baseurl <- "https://panelapp.agha.umccr.org/api/v1/panels/?format=json&page="

# find all pages through API call
panelapp_australia_last_page <- find_last_page(baseurl)

panelapp_australia_pages <- list()
for (i in 1:panelapp_australia_last_page){
  panelapp_australia_page <- fromJSON(paste0(baseurl, "&page=", i))
  message("Retrieving page ", i)
  panelapp_australia_pages[[i + 1]] <- panelapp_australia_page$results
}

# combine all into one
panelapp_australia <- rbind_pages(panelapp_australia_pages) %>%
  tibble() %>%
  select(id, name, version, version_created) %>%
  mutate(api_call = paste0("https://panelapp.agha.umccr.org/api/v1/panels/",
    id,
    "/?format=json")) %>%
  mutate(panel_source = "panelapp_australia")

# combine into one tibble
# filter for kidney related panels
filter_string <- config_vars_proj$panelap_filter_string

panelapp_panels <- bind_rows(panelapp_uk, panelapp_australia) %>%
  mutate(kidney_disease = str_detect(name, filter_string))

panelapp_panels_kidney <- panelapp_panels %>%
  filter(kidney_disease == TRUE)

############################################



############################################
## load all kidney related PanelApp panels into one tibble

# TODO: download API call data for reproducibility

panelapp_genes <- panelapp_panels_kidney %>%
  rowwise() %>%
  mutate(panel = list(fromJSON(api_call))) %>%
  unnest_wider(panel, names_repair = "unique") %>%
  select(id = id...1,
    name = name...2,
    version = version...3,
    version_created = version_created...4,
    panel_source, genes) %>%
  unnest(genes, names_repair = "unique") %>%
  unnest(gene_data, names_repair = "unique") %>%
  select(panel_id = id,
    panel_name = name,
    panel_version = version,
    panel_version_created = version_created,
    panel_source,
    entity_name,
    entity_type,
    gene_symbol,
    hgnc_id,
    omim_gene,
    evidence,
    confidence_level,
    mode_of_inheritance,
    phenotypes, publications) %>%
  ungroup()

# sapply used to collapse all list entries into one string
all_panelapp_genes <- panelapp_genes %>%
  arrange(panel_id, entity_name) %>%
  mutate(evidence = sapply(evidence, paste, collapse = "; "),
    phenotypes = sapply(phenotypes, paste, collapse = "; "),
    publications = sapply(publications, paste, collapse = "; ")) %>%
  mutate(evidence_category =
    case_when(
      confidence_level == 3 ~ "Green",
      confidence_level == 2 ~ "Amber",
      confidence_level == 1 | confidence_level == 0 ~ "Red",
    )
  ) %>%
  mutate(panel = paste0(panel_name,
      " (id=",
      panel_id,
      ", version=",
      panel_version,
      ", confidence=",
      evidence_category,
      ")")) %>%
  group_by(entity_name) %>%
  summarise(entity_name = paste(unique(entity_name), collapse = " | "),
    entity_type = paste(unique(entity_type), collapse = " | "),
    gene_symbol = paste(unique(gene_symbol), collapse = " | "),
    hgnc_id = paste(unique(hgnc_id), collapse = " | "),
    omim_gene = paste(unique(omim_gene), collapse = " | "),
    evidence = paste(unique(evidence), collapse = " | "),
    evidence_category = paste(unique(evidence_category), collapse = " | "),
    confidence_level = paste(unique(confidence_level), collapse = " | "),
    mode_of_inheritance = paste(unique(mode_of_inheritance), collapse = " | "),
    phenotypes = paste(unique(phenotypes), collapse = " | "),
    publications = paste(unique(publications), collapse = " | "),
    panel = paste(unique(panel), collapse = " | "),
    panel_id = paste(unique(panel_id), collapse = " | "),
    panel_name = paste(unique(panel_name), collapse = " | "),
    panel_version = paste(unique(panel_version), collapse = " | "),
    panel_version_created = paste(unique(panel_version_created),
      collapse = " | "),
    panel_source = paste(unique(panel_source), collapse = " | "),
    PanelApp_Green_or_amber = grepl("Green|Amber", evidence_category),
    source_count = n(),
    .groups = "keep") %>%
  ungroup() %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(gene_symbol)) %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id))

all_panelapp_genes_format <- all_panelapp_genes %>%
  select(approved_symbol,
    hgnc_id,
    gene_name_reported = gene_symbol,
    source = panel,
    source_count,
    source_evidence = PanelApp_Green_or_amber)

# TODO: normalize source_evidence to 0/1 as percentiles
# TODO: write a function for this normalization step

############################################


############################################
## save results
creation_date <- strftime(as.POSIXlt(Sys.time(),
  "UTC",
  "%Y-%m-%dT%H:%M:%S"), "%Y-%m-%d")

write_csv(panelapp_panels_kidney,
  file = paste0("results/01_PanelApp_panels.",
    creation_date,
    ".csv"),
  na = "NULL")

gzip(paste0("results/01_PanelApp_panels.", creation_date, ".csv"),
  overwrite = TRUE)

write_csv(all_panelapp_genes_format,
  file = paste0("results/01_PanelApp_genes.",
    creation_date,
    ".csv"),
  na = "NULL")

gzip(paste0("results/01_PanelApp_genes.", creation_date, ".csv"),
  overwrite = TRUE)
############################################
