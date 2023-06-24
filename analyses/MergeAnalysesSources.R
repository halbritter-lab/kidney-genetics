############################################
## load libraries
library(tidyverse)  ## needed for general table operations
library(readr)  ## needed to read files
library(tools)  ## needed for checksums
library("R.utils")  ## gzip downloaded and result files
library(config)
############################################


############################################
## define relative script path
project_topic <- "nephrology"
project_name <- "kidney-genetics"
script_path <- "/analyses/"

## read configs
config_vars_proj <- config::get(file = Sys.getenv("CONFIG_FILE"),
    config = project_topic)

## set working directory
setwd(paste0(config_vars_proj$projectsdir, project_name, script_path))

## set global options
options(scipen = 999)
############################################


############################################
## load all analyses files and transform table

# define analyses paths
analyses_paths <- c("01_PanelApp/results/",
  "02_Literature/results/",
  "03_DiagnosticPanels/results/",
  "04_HPO/results/",
  "05_PubTator/results/")

# find all CSV files in results folders and filter
# select only genes files
# select newest file
results_csv_table <- list.files(path = analyses_paths,
    pattern = ".csv",
    full.names = TRUE) %>%
  as_tibble() %>%
  separate(value, c("analysis", "path", "file"), sep = "\\/") %>%
  mutate(file_path = paste0(analysis, "/", path, "/", file)) %>%
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
  select(approved_symbol, hgnc_id, analysis, source_evidence) %>%
  group_by(approved_symbol, hgnc_id) %>%
  mutate(evidence_count = sum(source_evidence == TRUE)) %>%
  mutate(list_count =
    sum(source_evidence == TRUE | source_evidence == FALSE)) %>%
  ungroup %>%
  pivot_wider(
    names_from = analysis,
    values_from = source_evidence
  )

#TODO: annotate with OMIM P numbers
#TODO: annotate with GeneCC

#TODO: annotate with HPO kidney groups (cystic, nephrotic, cancer,...)
# annotate kidney groups
# https://clinicalgenome.org/working-groups/clinical-domain/clingen-kidney-disease-clinical-domain-working-group/
# 1) Complement-Mediated Kidney Diseases Gene Curation Expert Panel: https://www.clinicalgenome.org/affiliation/40069/ (no gene list, use https://panelapp.agha.umccr.org/panels/224/)
# 2) Congenital Anomalies of the Kidney and Urinary Tract Gene Curation Expert Panel: https://www.clinicalgenome.org/affiliation/40070/ (no gene list, use https://panelapp.genomicsengland.co.uk/panels/234/)
# 3) Glomerulopathy Gene Curation Expert Panel: https://www.clinicalgenome.org/affiliation/40068/
# 4) Kidney Cystic and Ciliopathy Disorders Gene Curation Expert Panel: https://www.clinicalgenome.org/affiliation/40066/
# 5) Tubulopathy Gene Curation Expert Panel: https://www.clinicalgenome.org/affiliation/40067/
# 6) Hereditary Cancer Gene Curation Expert Panel: https://clingen.info/affiliation/40023/

# general workflow:
# 1) get all genes from all kidney groups
# 2) get phenotype annotation table from HPO, group by gene and filter for kidney phenotypes
# 3) compute for each list in 1) the relative frequency of kidney phenotypes in 2), this nwill be the kidney group score
# 4) compute for each gene the a list of all kidney group scores and sort by the highest score (the gene is in the kidney group with the highest score, but this needs to be checked manually)
# 5) genes with score 0 are not in any kidney group and are annotated as "none"
# 6) genes with no annotated OMIM phenotypes are annotated as "NA"

# syndromic vs non-syndromic (categories in OMIM: GROWTH, SKELETAL, NEUROLOGIC, HEAD & NECK; exclude: CARDIOVASCULAR, ABDOMEN, GENITOURINARY)
# pediatric vs adult onset (OMIM: HPO terms Adult onset HP:0003581, Pediatric onset HP:0410280, maybe othe terms children of terms)

############################################


############################################
## save results
creation_date <- strftime(as.POSIXlt(Sys.time(),
  "UTC",
  "%Y-%m-%dT%H:%M:%S"), "%Y-%m-%d")

write_csv(results_genes_wider,
  file = paste0("merged/KidneyGenetics_MergeAnalysesSources.",
    creation_date,
    ".csv"),
  na = "NULL")

gzip(paste0("results/KidneyGenetics_MergeAnalysesSources.", creation_date, ".csv"),
  overwrite = TRUE)
############################################