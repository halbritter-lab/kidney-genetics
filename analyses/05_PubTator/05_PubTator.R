############################################
## load libraries
library(readr)        ## needed to read files
library(tidyverse)    ## needed for general table operations
library(httr)         ## needed for scraping
library(jsonlite)     ## needed for HGNC requests
library("R.utils")    ## gzip downloaded and result files
library(config)       ## needed for config loading
############################################


############################################
## define relative script path
project_topic <- "nephrology"
project_name <- "kidney-genetics"
script_path <- "/analyses/05_PubTator/"

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
## define API key
NCBI_API_KEY <- config_vars_proj$ncbi_api_key
############################################


############################################
# load global functions
# hgnc functions
source("../functions/hgnc-functions.R", local = TRUE)
# Pubtator functions
source("../functions/PubTator-functions.R", local = TRUE)
# NCBI functions
source("../functions/NCBI-datasets-v2-API-functions.R", local = TRUE)
# helper functions
source("../functions/helper-functions.R", local = TRUE)
############################################


############################################
## perform analysis

# define search query
# TODO: this query should be in a config file
search_query <- '("kidney disease" OR "renal disease") AND (gene OR syndrome) AND (variant OR mutation)'

# get number of pages for search query
pages_request <- pubtator_pages_request(search_query)

# make a tibble of all pages and perform request using the page as input
page_tibble <- 1:pages_request %>%
  as_tibble() %>%
  select(page = value) %>%
  rowwise() %>%
  mutate(results = list(pubtator_genes_in_request(search_query, page))) %>%
  ungroup()

page_tibble_unnest <- page_tibble %>%
  unnest(results) %>%
  filter(!is.na(text_identifier)) %>%
  separate_rows(text_identifier, sep = ";")

page_tibble_genes <- page_tibble_unnest %>%
  mutate(info = gene_info_from_gene_id(text_identifier)) %>%
  unnest(info)

# filter for human only, select columns
pubtator_genes <- page_tibble_genes %>%
  filter(tax_id == config_vars_proj$pubtator_tax_id) %>%
  select(pmid, symbol, text_part) %>%
  group_by(symbol) %>%
  summarise(symbol = paste(unique(symbol), collapse = "; "),
    pmid = paste(unique(paste0("PMID_", pmid)), collapse = " | "),
    text_part = paste(text_part, collapse = "; "),
    source_count = n(),
    .groups = "keep") %>%
  ungroup()

# format using same functions as in other analyses
pubtator_genes_format <- pubtator_genes %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(symbol)) %>%
  filter(!is.na(hgnc_id)) %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  mutate(at_least_three_pub = (source_count > 2)) %>%
  select(approved_symbol,
    hgnc_id,
    gene_name_reported = symbol,
    source = pmid,
    source_count,
    source_evidence = at_least_three_pub)

# normalize source_count to 0/1 as percentiles
pubtator_genes_format_normalize <- pubtator_genes_format %>%
  normalize_percentile("source_count")
############################################


############################################
## save results
write_csv(pubtator_genes_format_normalize,
  file = paste0("results/05_PubTator_genes.",
    current_date,
    ".csv"),
  na = "NULL")

gzip(paste0("results/05_PubTator_genes.", current_date, ".csv"),
  overwrite = TRUE)
############################################
