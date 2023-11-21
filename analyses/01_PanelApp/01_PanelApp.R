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

# compute date only once or somehow in config
current_date <- get_current_date_iso8601()
############################################


############################################
# load global functions
# hgnc functions
source("../functions/hgnc-functions.R", local = TRUE)
source("../functions/panelapp-functions.R", local = TRUE)
source("../functions/file-functions.R", local = TRUE)
# helper functions
source("../functions/helper-functions.R", local = TRUE)
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
filter_string <- config_vars_proj$panelapp_filter_string

panelapp_panels <- bind_rows(panelapp_uk, panelapp_australia) %>%
  mutate(kidney_disease = str_detect(name, filter_string))

panelapp_panels_kidney <- panelapp_panels %>%
  filter(kidney_disease == TRUE)

############################################


############################################
## load all kidney related PanelApp panels into one tibble

# directory to save JSON files
download_path <- "data/downloads"

# check if files are older than 1 month
panelapp_panels_kidney_checked <- panelapp_panels_kidney %>%
  mutate(json_file_name = paste0(id, "_", str_replace_all(name, "[[:punct:]]", ""))) %>%
  rowwise() %>%
  mutate(file_is_young = check_file_age(json_file_name, download_path, 1))

# panelapp JSON file download
if (all(panelapp_panels_kidney_checked$file_is_young)) {
  panelapp_panels_kidney_files <- panelapp_panels_kidney_checked %>%
    rowwise() %>%
    mutate(json_download_file = get_newest_file(json_file_name, download_path))
} else {
  panelapp_panels_kidney_files <- panelapp_panels_kidney_checked %>%
    mutate(json_file_path = file.path(download_path, paste0(json_file_name, ".json"))) %>%
    rowwise() %>%
    mutate(json_download_file = download_and_save_json(api_call, json_file_path))
}

# now read the local JSON files using jsonlite and reprocess
panelapp_genes <- panelapp_panels_kidney_files %>%
  rowwise() %>%
  mutate(panel = list(fromJSON(json_download_file))) %>%
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
  mutate(
    evidence = sapply(evidence, paste, collapse = "; "),
    phenotypes = sapply(phenotypes, paste, collapse = "; "),
    publications = sapply(publications, paste, collapse = "; ")
  ) %>%
  mutate(evidence_category =
    case_when(
      confidence_level == 3 ~ "Green",
      confidence_level == 2 ~ "Amber",
      confidence_level == 1 | confidence_level == 0 ~ "Red"
    )
  ) %>%
  mutate(
    panel = paste0(panel_name,
      " (id=",
      panel_id,
      ", version=",
      panel_version,
      ", confidence=",
      evidence_category,
      ")"),
    Green_or_Amber_count = as.integer(grepl("Green|Amber", evidence_category))
  ) %>%
  group_by(entity_name) %>%
  summarise(
    entity_name = paste(unique(entity_name), collapse = " | "),
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
    panel_version_created = paste(unique(panel_version_created), collapse = " | "),
    panel_source = paste(unique(panel_source), collapse = " | "),
    source_count = sum(Green_or_Amber_count),
    .groups = "keep"
  ) %>%
  ungroup() %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(gene_symbol)) %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  mutate(source_evidence = grepl("Green|Amber", evidence_category))

all_panelapp_genes_format <- all_panelapp_genes %>%
  select(
    approved_symbol,
    hgnc_id,
    gene_name_reported = gene_symbol,
    source = panel,
    source_count,
    source_evidence
  )

# normalize source_count to 0/1 as percentiles
all_panelapp_genes_format_normalize <- all_panelapp_genes_format %>%
  normalize_percentile("source_count")

############################################


############################################
## save results
write_csv(panelapp_panels_kidney,
  file = paste0("results/01_PanelApp_panels.",
    current_date,
    ".csv"),
  na = "NULL")

gzip(paste0("results/01_PanelApp_panels.", current_date, ".csv"),
  overwrite = TRUE)

write_csv(all_panelapp_genes_format_normalize,
  file = paste0("results/01_PanelApp_genes.",
    current_date,
    ".csv"),
  na = "NULL")

gzip(paste0("results/01_PanelApp_genes.", current_date, ".csv"),
  overwrite = TRUE)
############################################
