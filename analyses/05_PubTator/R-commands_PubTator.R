############################################
## load libraries
library(readr)
library(tidyverse)
library(httr)
library(jsonlite)
library(config)
############################################


############################################
## define relativescript path
project_name <- "kidney-genetics"
script_path <- "/analyses/05_PubTator/"
## read config
config_vars <- config::get(file = Sys.getenv("CONFIG_FILE"))
## set working directory
setwd(paste0(config_vars$projectsdir, project_name, script_path))
## set global options
options(scipen = 999)
############################################


############################################
## define API key
# TODO: this has to go into config
NCBI_API_KEY <- config_vars$ncbi_api_key
############################################


############################################
# load global functions
# hgnc functions 
source("../functions/hgnc-functions.R", local = TRUE)
# NCBI functions 
source("../functions/NCBI-datasets-v2-API-functions.R", local = TRUE)
############################################


############################################
## perform analysis

# define search query
search_query <- '("kidney disease" OR "renal disease") AND (gene OR syndrome) AND (variant OR mutation)'

# get number of pages for search query
pages_request <- pubtator_pages_request(search_query)

# make a tibble of all pages and perform request using thhe page as input
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
  filter(tax_id == 9606) %>%
  select(pmid, symbol, text_part) %>%
  group_by(symbol) %>%
  summarise(symbol = paste(unique(symbol), collapse = "; "),
    pmid = paste(unique(pmid), collapse = "; "),
    text_part = paste(text_part, collapse = "; "),
    source_count = n(),
    .groups = "keep") %>%
  ungroup()

# normalize using same functions as in other analyses
pubtator_genes_normalize <- pubtator_genes %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(symbol)) %>%
  filter(!is.na(hgnc_id)) %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  select(approved_symbol, hgnc_id, gene_name_reported = symbol, source = pmid, source_count, source_evidence = text_part)

# TODO: normalize source_evidence as in other analyses (maybe pub count > 2)

############################################


############################################
## save results
creation_date <- strftime(as.POSIXlt(Sys.time(), "UTC", "%Y-%m-%dT%H:%M:%S"), "%Y-%m-%d")
write_csv(pubtator_genes_normalize, file = paste0("results/05_PubTator_genes.", creation_date,".csv"), na = "NULL")
############################################
