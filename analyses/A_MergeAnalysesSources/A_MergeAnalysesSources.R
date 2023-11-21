############################################
## load libraries
library(tidyverse)  ## needed for general table operations
library(readr)      ## needed to read files
library(tools)      ## needed for checksums
library("R.utils")  ## gzip downloaded and result files
library(config)     ## needed for config loading
############################################


############################################
## define relative script path
project_topic <- "nephrology"
project_name <- "kidney-genetics"
script_path <- "/analyses/A_MergeAnalysesSources/"

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
# helper functions
source("../functions/helper-functions.R", local = TRUE)
############################################


############################################
# compute date only once or somehow in config
current_date <- get_current_date_iso8601()
############################################


############################################
## load all analyses files and transform table

# define analyses paths
analyses_paths <- c("/analyses/01_PanelApp/results/",
  "/analyses/02_Literature/results/",
  "/analyses/03_DiagnosticPanels/results/",
  "/analyses/04_HPO/results/",
  "/analyses/05_PubTator/results/")

# find all CSV files in results folders and filter
# select only genes files
# select newest file
results_csv_table <- list.files(path = paste0(config_vars_proj$projectsdir, project_name, analyses_paths),
    pattern = ".csv.gz",
    full.names = TRUE) %>%
  as_tibble() %>%
  mutate(file_path = str_replace_all(value, paste0(config_vars_proj$projectsdir, project_name, "/analyses/"), "\\.\\.\\/")) %>%
  mutate(value = str_replace_all(value, paste0(config_vars_proj$projectsdir, project_name, "/analyses/"), "")) %>%
  separate(value, c("analysis", "path", "file"), sep = "\\/") %>%
  separate(file, c(NA, "analysis_date", NA), sep = "\\.") %>%
  mutate(results_file_id = row_number()) %>%
  mutate(md5sum_file = md5sum(file_path)) %>%
  dplyr::select(results_file_id,
    file_path,
    analysis,
    analysis_date,
    md5sum_file) %>%
  filter(str_detect(file_path, "genes")) %>%
  group_by(analysis) %>%
    filter(analysis_date == max(analysis_date)) %>%
  ungroup() %>%
  arrange(analysis)

# load the csv files
results_genes <- results_csv_table %>%
  rowwise() %>%
  mutate(genes_list = list(read_csv(file_path,
    na = "NULL",
    col_types = cols(approved_symbol = col_character(),
      hgnc_id = col_character(), gene_name_reported = col_character(),
      source = col_character(), source_count = col_double(),
      source_evidence = col_logical())
    ))) %>%
  ungroup() %>%
  select(analysis, genes_list) %>%
  unnest(genes_list)

# generate wide table and compute
# evidence_count = sum of lists where the source_evidence is TRUE
# list_count = sum lists where gene is found (source_evidence is TRUE or FALSE)
results_genes_wider <- results_genes %>%
  select(approved_symbol, hgnc_id, analysis, source_evidence, source_count_percentile) %>%
  group_by(approved_symbol, hgnc_id) %>%
  mutate(evidence_count = sum(source_evidence == TRUE),
    source_count_percentile = sum(source_count_percentile)) %>%
  mutate(list_count =
    sum(source_evidence == TRUE | source_evidence == FALSE)) %>%
  ungroup %>%
  pivot_wider(
    names_from = analysis,
    values_from = source_evidence
  )

############################################


############################################
## save results
write_csv(results_genes_wider,
  file = paste0("results/A_MergeAnalysesSources.",
    current_date,
    ".csv"),
  na = "NULL")

gzip(paste0("results/A_MergeAnalysesSources.", current_date, ".csv"),
  overwrite = TRUE)

write_csv(results_csv_table,
  file = paste0("results/A_results_csv_table.",
    current_date,
    ".csv"),
  na = "NULL")

gzip(paste0("results/A_results_csv_table.", current_date, ".csv"),
  overwrite = TRUE)
############################################