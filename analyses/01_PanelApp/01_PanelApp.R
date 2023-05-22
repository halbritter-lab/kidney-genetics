############################################
## load libraries
library(readr)
library(tidyverse)
library(rvest)
library(jsonlite)
library(config)
############################################


############################################
## define relative script path
project_name <- "kidney-genetics"
script_path <- "/analyses/01_PanelApp/"
## read config
config_vars <- config::get(file = Sys.getenv("CONFIG_FILE"))
## set working directory
setwd(paste0(config_vars$projectsdir, project_name, script_path))
## set global options
options(scipen = 999)
############################################


############################################
# load global functions
# hgnc functions
source("../functions/hgnc-functions.R", local = TRUE)
############################################


############################################
## get a list of all PanelApp panels and filter

## PanelApp UK
# store all pages in a list first
baseurl <- "https://panelapp.genomicsengland.co.uk/api/v1/panels/?format=json"
panelapp_uk_pages <- list()
for (i in 1:4){
  panelapp_uk_page <- fromJSON(paste0(baseurl, "&page=", i))
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
  mutate(panel_source = "panelapp_uk") %>%
  filter(str_detect(name, "[Kk]idney|[Rr]enal|[Nn]ephro"))

## PanelApp Australia
# store all pages in a list first
baseurl <- "https://panelapp.agha.umccr.org/api/v1/panels/?format=json"
panelapp_australia_pages <- list()
for (i in 1:3){
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
  mutate(panel_source = "panelapp_australia") %>%
  filter(str_detect(name, "[Kk]idney|[Rr]enal|[Nn]ephro"))

# combine into one tibble
panelapp_panels <- bind_rows(panelapp_uk, panelapp_australia)

############################################



############################################
## load all kidney related PanelApp panels into one tibble

# TODO: download API call data for reproducibility

panelapp_genes <- panelapp_panels %>%
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

all_panelapp_genes <- panelapp_genes %>%
  mutate(evidence = sapply(evidence, paste, collapse = "; "),
    phenotypes = sapply(phenotypes, paste, collapse = "; "),
    publications = sapply(publications, paste, collapse = "; ")) %>%
  group_by(entity_name) %>%
  summarise(entity_name = paste(unique(entity_name), collapse = " | "),
    entity_type = paste(unique(entity_type), collapse = " | "),
    gene_symbol = paste(unique(gene_symbol), collapse = " | "),
    hgnc_id = paste(unique(hgnc_id), collapse = " | "),
    omim_gene = paste(unique(omim_gene), collapse = " | "),
    evidence = paste(unique(evidence), collapse = " | "),
    confidence_level = paste(unique(confidence_level), collapse = " | "),
    mode_of_inheritance = paste(unique(mode_of_inheritance), collapse = " | "),
    phenotypes = paste(unique(phenotypes), collapse = " | "),
    publications = paste(unique(publications), collapse = " | "),
    panel = paste(paste0(panel_name,
      " (id=",
      panel_id,
      ", version=",
      panel_version,
      ")"), collapse = " | "),
    panel_id = paste(unique(panel_id), collapse = " | "),
    panel_name = paste(unique(panel_name), collapse = " | "),
    panel_version = paste(unique(panel_version), collapse = " | "),
    panel_version_created = paste(unique(panel_version_created),
      collapse = " | "),
    panel_source = paste(unique(panel_source), collapse = " | "),
    PanelApp_Green_or_amber = grepl("Green|Amber", evidence),
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
# TODO: add green/amber/red status to source

############################################


############################################
## save results
# TODO: gzip csv result files
creation_date <- strftime(as.POSIXlt(Sys.time(),
  "UTC",
  "%Y-%m-%dT%H:%M:%S"), "%Y-%m-%d")

write_csv(panelapp_panels,
  file = paste0("results/01_PanelApp_panels.",
    creation_date,
    ".csv"),
  na = "NULL")

write_csv(all_panelapp_genes_format,
  file = paste0("results/01_PanelApp_genes.",
    creation_date,
    ".csv"),
  na = "NULL")
############################################
