############################################
## load libraries
library(tidyverse)  ## needed for general table operations
library(jsonlite)  ##needed for HGNC requests
library(rvest)  ## needed for scraping
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
# load global functions
# hgnc functions
source("functions/hgnc-functions.R", local = TRUE)
source("functions/hpo-functions.R", local = TRUE)
source("functions/file-functions.R", local = TRUE)
source("functions/gtex-functions.R", local = TRUE)
############################################


############################################
## load all analyses files and transform table

# define analyses paths
merged_path <- "merged/"

# find all CSV files in merged folders and filter
# select newest file
# TODO: replace with file function
merged_csv_table <- list.files(path = merged_path,
    pattern = ".csv",
    full.names = TRUE) %>%
  as_tibble() %>%
  separate(value, c("path", "file"), sep = "\\/") %>%
  mutate(file_path = paste0(path, "/", file)) %>%
  separate(file, c(NA, "file_date", NA), sep = "\\.") %>%
  mutate(results_file_id = row_number()) %>%
  mutate(md5sum_file = md5sum(file_path)) %>%
  dplyr::select(results_file_id,
    file_path,
    file_date,
    md5sum_file) %>%
  filter(file_date == max(file_date))

# load the csv files
merged_genes <- merged_csv_table %>%
  rowwise() %>%
  mutate(merged_list = list(read_csv(file_path,
    na = "NULL"
    ))) %>%
  ungroup() %>%
  select(merged_list) %>%
  unnest(merged_list)
############################################


############################################
## get relevant HPO terms for kidney disease classification
# get the current date
# TODO: compute date only once or somehow in config
current_date <- strftime(as.POSIXlt(Sys.time(), "UTC", "%Y-%m-%dT%H:%M:%S"), "%Y-%m-%d")

# 1) get all children of term upper urinary tract (HP:0010935) for classification into kidney disease groups
# walk through the ontology tree and add all unique terms descending from
# Abnormality of the upper urinary tract (HP:0010935)
# we load and use the results of previous walks through the ontology tree if not older then 1 month

if (check_file_age("hpo_list_kidney", "shared/", 1)) {
  hpo_list_kidney <- read_csv(get_newest_file("hpo_list_kidney", "shared"))
} else {
  all_hpo_children_list_kidney <- HPO_all_children_from_term("HP:0010935")

  # transform the list into a tibble
  hpo_list_kidney <- all_hpo_children_list_kidney %>%
    unlist() %>%
    tibble(`term` = .) %>%
    unique() %>%
    mutate(query_date = current_date)

  write_csv(hpo_list_kidney,
    file = paste0("shared/hpo_list_kidney.",
      current_date,
      ".csv"),
    na = "NULL")

  gzip(paste0("shared/hpo_list_kidney.", current_date, ".csv"),
    overwrite = TRUE)
}

# 2) get all children of term Adult onset (HP:0003581) for classification into adult vs pediatric onset
# we load and use the results of previous walks through the ontology tree if not older then 1 month

if (check_file_age("hpo_list_adult", "shared/", 1)) {
  hpo_list_adult <- read_csv(get_newest_file("hpo_list_adult", "shared"))
  hpo_list_non_adult <- read_csv(get_newest_file("hpo_list_non_adult", "shared"))
} else {
  # walk through the ontology tree and add all unique terms descending from
  # Adult onset (HP:0003581)
  all_hpo_children_list_adult <- HPO_all_children_from_term("HP:0003581")

  # transform the list into a tibble
  hpo_list_adult <- all_hpo_children_list_adult %>%
    unlist() %>%
    tibble(`term` = .) %>%
    unique() %>%
    mutate(query_date = current_date)

    write_csv(hpo_list_adult,
      file = paste0("shared/hpo_list_adult.",
        current_date,
        ".csv"),
      na = "NULL")

    gzip(paste0("shared/hpo_list_adult.", current_date, ".csv"),
      overwrite = TRUE)

  # walk through the ontology tree and add all unique terms descending from
  # Onset HP:0003674
  all_hpo_children_list_onset <- HPO_all_children_from_term("HP:0003674")

  # then remove all terms that are children of Adult onset (HP:0003581) to get
  # all terms that are children of Pediatric onset (HP:0410280), etc.
  # and transform the list into a tibble
  hpo_list_non_adult <- setdiff(all_hpo_children_list_onset,
      all_hpo_children_list_adult) %>%
    unlist() %>%
    tibble(`term` = .) %>%
    unique() %>%
    mutate(query_date = current_date)

    write_csv(hpo_list_non_adult,
      file = paste0("shared/hpo_list_non_adult.",
        current_date,
        ".csv"),
      na = "NULL")

    gzip(paste0("shared/hpo_list_non_adult.", current_date, ".csv"),
      overwrite = TRUE)
}

# 3) syndromic vs non-syndromic (categories in OMIM: GROWTH, SKELETAL, NEUROLOGIC, HEAD & NECK; exclude: CARDIOVASCULAR, ABDOMEN, GENITOURINARY)
# we load and use the results of previous walks through the ontology tree if not older then 1 month

if (check_file_age("hpo_list_syndromic", "shared/", 1)) {
  # load the list of syndromic terms
  hpo_list_syndromic <- read_csv(get_newest_file("hpo_list_syndromic", "shared"))
  # also load the sublists if not older then 1 month
  all_hpo_children_list_growth_term <- read_csv(get_newest_file("hpo_list_growth_term", "shared"))
  all_hpo_children_list_skeletal_term <- read_csv(get_newest_file("hpo_list_skeletal_term", "shared"))
  all_hpo_children_list_neurologic_term <- read_csv(get_newest_file("hpo_list_neurologic_term", "shared"))
  all_hpo_children_list_head_and_neck_term <- read_csv(get_newest_file("hpo_list_head_and_neck_term", "shared"))
} else {
  # walk through the ontology tree and add all unique terms descending from
  # Growth abnormality (HP:0001507)
  all_hpo_children_list_growth <- HPO_all_children_from_term("HP:0001507")

  # walk through the ontology tree and add all unique terms descending from
  # Skeletal system abnormality (HP:0000924)
  all_hpo_children_list_skeletal <- HPO_all_children_from_term("HP:0000924")

  # walk through the ontology tree and add all unique terms descending from
  # Neurologic abnormality (HP:0000707)
  all_hpo_children_list_neurologic <- HPO_all_children_from_term("HP:0000707")

  # walk through the ontology tree and add all unique terms descending from
  # Head and neck abnormality (HP:0000152)
  all_hpo_children_list_head_and_neck <- HPO_all_children_from_term("HP:0000152")

  # add the base term name and term id to each list
  # safe the lists as csv and gzip them
  all_hpo_children_list_growth_term <- all_hpo_children_list_growth %>%
    mutate(base_term = "Growth abnormality",
      base_term_id = "HP:0001507") %>%
    mutate(query_date = current_date) %>%
    select(term = value, base_term, base_term_id)

  all_hpo_children_list_skeletal_term <- all_hpo_children_list_skeletal %>%
    mutate(base_term = "Skeletal system abnormality",
      base_term_id = "HP:0000924") %>%
    mutate(query_date = current_date) %>%
    select(term = value, base_term, base_term_id)

  all_hpo_children_list_neurologic_term <- all_hpo_children_list_neurologic %>%
    mutate(base_term = "Neurologic abnormality",
      base_term_id = "HP:0000707") %>%
    mutate(query_date = current_date) %>%
    select(term = value, base_term, base_term_id)

  all_hpo_children_list_head_and_neck_term <- all_hpo_children_list_head_and_neck %>%
    mutate(base_term = "Head and neck abnormality",
      base_term_id = "HP:0000152") %>%
    mutate(query_date = current_date) %>%
    select(term = value, base_term, base_term_id)

  write_csv(all_hpo_children_list_growth_term,
    file = paste0("shared/hpo_list_growth_term.",
      current_date,
      ".csv"),
    na = "NULL")

  write_csv(all_hpo_children_list_skeletal_term,
    file = paste0("shared/hpo_list_skeletal_term.",
      current_date,
      ".csv"),
    na = "NULL")

  write_csv(all_hpo_children_list_neurologic_term,
    file = paste0("shared/hpo_list_neurologic_term.",
      current_date,
      ".csv"),
    na = "NULL")

  write_csv(all_hpo_children_list_head_and_neck_term,
    file = paste0("shared/hpo_list_head_and_neck_term.",
      current_date,
      ".csv"),
    na = "NULL")

  gzip(paste0("shared/hpo_list_growth_term.", current_date, ".csv"),
    overwrite = TRUE)

  gzip(paste0("shared/hpo_list_skeletal_term.", current_date, ".csv"),
    overwrite = TRUE)

  gzip(paste0("shared/hpo_list_neurologic_term.", current_date, ".csv"),
    overwrite = TRUE)

  gzip(paste0("shared/hpo_list_head_and_neck_term.", current_date, ".csv"),
    overwrite = TRUE)

  # merge all lists and transform the list into a tibble
  # summarize the list to unique terms
  # safe the list as csv and gzip it
  hpo_list_syndromic <- bind_rows(all_hpo_children_list_growth_term,
      all_hpo_children_list_skeletal_term,
      all_hpo_children_list_neurologic_term,
      all_hpo_children_list_head_and_neck_term) %>%
    group_by(term) %>%
    summarise(base_term = paste0(base_term, collapse = "; "),
      base_term_id = paste0(base_term_id, collapse = "; "),
      query_date = current_date)

  write_csv(hpo_list_syndromic,
    file = paste0("shared/hpo_list_syndromic.",
      current_date,
      ".csv"),
    na = "NULL")

  gzip(paste0("shared/hpo_list_syndromic.", current_date, ".csv"),
    overwrite = TRUE)
}
############################################


############################################
## download all required database sources from HPO and OMIM
# we load and use the results of previous walks through the ontology tree if not older then 1 month

# disease ontology annotations from HPO
if (check_file_age("phenotype", "shared/data/downloads/", 1)) {
  phenotype_hpoa_filename <- get_newest_file("phenotype", "shared/data/downloads/")
} else {
  # TODO: compute date only once or somehow in config
  file_date <- strftime(as.POSIXlt(Sys.time(), "UTC", "%Y-%m-%dT%H:%M:%S"), "%Y-%m-%d")

  # disease ontology annotations from HPO
  phenotype_hpoa_url <- config_vars_proj$phenotype_hpoa

  phenotype_hpoa_filename <- paste0("shared/data/downloads/phenotype.",
    file_date,
    ".hpoa")

  download.file(phenotype_hpoa_url, phenotype_hpoa_filename, mode = "wb")

  gzip(phenotype_hpoa_filename,
    overwrite = TRUE)

  phenotype_hpoa_filename <- paste0(phenotype_hpoa_filename,
    ".gz")
}

# OMIM links to genemap2 file needs to be set in config and applied for at
# https://www.omim.org/downloads
if (check_file_age("omim_genemap2", "shared/data/downloads/", 1)) {
  omim_genemap2_filename <- get_newest_file("omim_genemap2", "shared/data/downloads/")
} else {
  # TODO: compute date only once or somehow in config
  file_date <- strftime(as.POSIXlt(Sys.time(), "UTC", "%Y-%m-%dT%H:%M:%S"), "%Y-%m-%d")

  # OMIM links to genemap2 file needs to be set in config and applied for at
  # https://www.omim.org/downloads
  omim_genemap2_url <- config_vars_proj$omim_genemap2_url

  omim_genemap2_filename <- paste0("shared/data/downloads/omim_genemap2.",
    file_date,
    ".txt")

  download.file(omim_genemap2_url, omim_genemap2_filename, mode = "wb")

  gzip(omim_genemap2_filename,
    overwrite = TRUE)

  omim_genemap2_filename <- paste0(omim_genemap2_filename,
    ".txt")
}
############################################


############################################
## load OMIM and terms files and reformat them
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
# annotate kidney groups from:
# https://clinicalgenome.org/working-groups/clinical-domain/clingen-kidney-disease-clinical-domain-working-group/

# TODO: safe the results of the API calls or websites for reproducibility

## general workflow:
# A) get all genes from all kidney groups
# B) get phenotype annotation table from HPO, group by gene and filter for kidney phenotypes
# C) compute for each list in A) the relative frequency of kidney phenotypes in 2), this will be the kidney group score
# D) compute for each gene a list of all kidney group scores and sort by the highest score (the gene is in the kidney group with the highest score, but this needs to be checked manually)

## A get all genes from all kidney groups
# 1) Complement-Mediated Kidney Diseases Gene Curation Expert Panel: https://www.clinicalgenome.org/affiliation/40069/ (no gene list, use https://panelapp.agha.umccr.org/panels/224/)

api_call <- paste0(
    "https://panelapp.agha.umccr.org/api/v1/panels/",
    "224",
    "/?format=json")

complement_mediated_kidney_diseases_response <- fromJSON(api_call)

complement_mediated_kidney_diseases_gene_list <- complement_mediated_kidney_diseases_response$gene %>%
  tibble() %>%
  unnest_wider(gene_data, names_repair = "unique") %>%
  select(entity_name,
    entity_type,
    gene_symbol,
    hgnc_id,
    omim_gene,
    evidence,
    confidence_level,
    mode_of_inheritance,
    phenotypes, publications) %>%
  mutate(evidence = sapply(evidence, paste, collapse = "; "),
    phenotypes = sapply(phenotypes, paste, collapse = "; "),
    publications = sapply(publications, paste, collapse = "; ")) %>%
  filter(confidence_level == 3) %>%
  select(gene_symbol) %>%
  distinct() %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(gene_symbol)) %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  mutate(kidney_disease_group = "complement_mediated_kidney_diseases") %>%
  mutate(kidney_disease_group_short = "complement") %>%
  select(approved_symbol,
    hgnc_id,
    gene_name_reported = gene_symbol,
    kidney_disease_group,
    kidney_disease_group_short)

# 2) Congenital Anomalies of the Kidney and Urinary Tract Gene Curation Expert Panel: https://www.clinicalgenome.org/affiliation/40070/ (no gene list, use https://panelapp.genomicsengland.co.uk/panels/234/)

api_call <- paste0(
    "https://panelapp.genomicsengland.co.uk/api/v1/panels/",
    "234",
    "/?format=json")

congenital_anomalies_of_the_kidney_and_urinary_tract_response <- fromJSON(api_call)

congenital_anomalies_of_the_kidney_and_urinary_tract_gene_list <- congenital_anomalies_of_the_kidney_and_urinary_tract_response$gene %>%
  tibble() %>%
  unnest_wider(gene_data, names_repair = "unique") %>%
  select(entity_name,
    entity_type,
    gene_symbol,
    hgnc_id,
    omim_gene,
    evidence,
    confidence_level,
    mode_of_inheritance,
    phenotypes, publications) %>%
  mutate(evidence = sapply(evidence, paste, collapse = "; "),
    phenotypes = sapply(phenotypes, paste, collapse = "; "),
    publications = sapply(publications, paste, collapse = "; ")) %>%
  filter(confidence_level == 3) %>%
  select(gene_symbol) %>%
  distinct() %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(gene_symbol)) %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  mutate(kidney_disease_group = "congenital_anomalies_of_the_kidney_and_urinary_tract") %>%
  mutate(kidney_disease_group_short = "cakut") %>%
  select(approved_symbol,
    hgnc_id,
    gene_name_reported = gene_symbol,
    kidney_disease_group,
    kidney_disease_group_short)

# 3) Glomerulopathy Gene Curation Expert Panel: https://www.clinicalgenome.org/affiliation/40068/

url <- "https://www.clinicalgenome.org/affiliation/40068/"

glomerulopathy_page <- read_html(url)

glomerulopathy_gene_list <- glomerulopathy_page %>%
    html_nodes(xpath = '//p[contains(text(),"Gene list:")]') %>%
  html_text() %>%
  str_remove_all("Gene list: ") %>%
  tibble(`gene_symbol` = .) %>%
  separate_rows(., gene_symbol, convert = TRUE) %>%
  select(gene_symbol) %>%
  distinct() %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(gene_symbol)) %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  mutate(kidney_disease_group = "glomerulopathy") %>%
  mutate(kidney_disease_group_short = "glomerulopathy") %>%
  select(approved_symbol,
    hgnc_id,
    gene_name_reported = gene_symbol,
    kidney_disease_group,
    kidney_disease_group_short)

# 4) Kidney Cystic and Ciliopathy Disorders Gene Curation Expert Panel: https://www.clinicalgenome.org/affiliation/40066/

url <- "https://www.clinicalgenome.org/affiliation/40066/"

kidney_cystic_and_ciliopathy_disorders_page <- read_html(url)

kidney_cystic_and_ciliopathy_disorders_gene_list <- kidney_cystic_and_ciliopathy_disorders_page %>%
    html_nodes(xpath = '//em[contains(text(),"ADAMTS9")]') %>%
  html_text() %>%
  tibble(`gene_symbol` = .) %>%
  separate_rows(., gene_symbol, convert = TRUE) %>%
  select(gene_symbol) %>%
  distinct() %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(gene_symbol)) %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  mutate(kidney_disease_group = "kidney_cystic_and_ciliopathy_disorders") %>%
  mutate(kidney_disease_group_short = "cyst_cilio") %>%
  select(approved_symbol,
    hgnc_id,
    gene_name_reported = gene_symbol,
    kidney_disease_group,
    kidney_disease_group_short)

# 5) Tubulopathy Gene Curation Expert Panel: https://www.clinicalgenome.org/affiliation/40067/

url <- "https://www.clinicalgenome.org/affiliation/40067/"

tubulopathy_page <- read_html(url)

tubulopathy_gene_list <- tubulopathy_page %>%
    html_nodes(xpath = '//p[contains(text(),"Gene list:")]') %>%
  html_text() %>%
  str_remove_all("Gene list: ") %>%
  str_replace_all("CLCNKB", "CLCNKB") %>% #fix for typo in html
  tibble(`gene_symbol` = .) %>%
  separate_rows(., gene_symbol, convert = TRUE) %>%
  select(gene_symbol) %>%
  distinct() %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(gene_symbol)) %>%
  filter(!is.na(hgnc_id)) %>% #fix for typo in html
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  mutate(kidney_disease_group = "tubulopathy") %>%
  mutate(kidney_disease_group_short = "tubulopathy") %>%
  select(approved_symbol,
    hgnc_id,
    gene_name_reported = gene_symbol,
    kidney_disease_group,
    kidney_disease_group_short)

# 6) Hereditary Cancer Gene Curation Expert Panel: https://clingen.info/affiliation/40023/
# API call: https://search.clinicalgenome.org/api/affiliates/10023?queryParams

api_call <- "https://search.clinicalgenome.org/api/affiliates/10023?queryParams"

hereditary_cancer_response <- fromJSON(api_call)

hereditary_cancer_gene_list <- hereditary_cancer_response$rows %>%
  tibble() %>%
  select(gene_symbol = symbol) %>%
  distinct() %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(gene_symbol)) %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  mutate(kidney_disease_group = "hereditary_cancer") %>%
  mutate(kidney_disease_group_short = "cancer") %>%
  select(approved_symbol,
    hgnc_id,
    gene_name_reported = gene_symbol,
    kidney_disease_group,
    kidney_disease_group_short)

## bind all tables
all_kidney_groups <- bind_rows(complement_mediated_kidney_diseases_gene_list,
    congenital_anomalies_of_the_kidney_and_urinary_tract_gene_list,
    glomerulopathy_gene_list,
    kidney_cystic_and_ciliopathy_disorders_gene_list,
    tubulopathy_gene_list,
    hereditary_cancer_gene_list)

# B) get phenotype annotation table from HPO, group by gene and filter for kidney phenotypes
# C) compute for each list in A) the relative frequency of kidney phenotypes in 2), this will be the kidney group score
# D) compute for each gene a list of all kidney group scores and sort by the highest score (the gene is in the kidney group with the highest score, but this needs to be checked manually)
# TODO: describe the logic of the analysis in comments step by step
omim_genemap2_disease_and_gene <- omim_genemap2 %>%
  select(disease_ontology_id, approved_symbol)

phenotype_hpoa_filter_kidney <- phenotype_hpoa %>%
   filter(hpo_id %in% hpo_list_kidney$term) %>%
   select(database_id, hpo_id) %>%
   unique()

hpo_gene_list_kidney <- phenotype_hpoa_filter_kidney %>%
  left_join(omim_genemap2_disease_and_gene, by = c("database_id" = "disease_ontology_id"), relationship = "many-to-many") %>%
  filter(!is.na(approved_symbol)) %>%
  unique()

all_kidney_groups_hpo <- all_kidney_groups %>%
  left_join(hpo_gene_list_kidney, by = c("approved_symbol"), relationship = "many-to-many") %>%
  filter(!is.na(hpo_id)) %>%
  select(kidney_disease_group_short, hpo_id) %>%
  group_by(kidney_disease_group_short) %>%
  mutate(hpo_id_count = n()) %>%
  ungroup() %>%
  group_by(kidney_disease_group_short, hpo_id) %>%
  summarise(hpo_id = unique(hpo_id),
    hpo_id_count = unique(hpo_id_count),
    hpo_id_group_count = n(),
    .groups = "keep") %>%
  ungroup() %>%
  mutate(hpo_id_group_frac = hpo_id_group_count / hpo_id_count) %>%
  select(-hpo_id_count, -hpo_id_group_count)

hpo_gene_list_kidney_all_kidney_groups <- hpo_gene_list_kidney %>%
  left_join(all_kidney_groups_hpo, by = c("hpo_id"), relationship = "many-to-many") %>%
  filter(!is.na(kidney_disease_group_short)) %>%
  select(approved_symbol, database_id , kidney_disease_group_short, hpo_id_group_frac) %>%
  unique() %>%
  group_by(approved_symbol, database_id, kidney_disease_group_short) %>%
  summarise(approved_symbol = unique(approved_symbol),
    database_id = unique(database_id),
    kidney_disease_group_short = unique(kidney_disease_group_short),
    hpo_id_group_p = round(sum(hpo_id_group_frac), 3),
    .groups = "keep") %>%
  ungroup() %>%
  arrange(approved_symbol, desc(hpo_id_group_p))

hpo_gene_list_kidney_all_kidney_groups_summarized_for_join <- hpo_gene_list_kidney_all_kidney_groups %>%
  group_by(approved_symbol, database_id) %>%
  summarise(approved_symbol = unique(approved_symbol),
    database_id = unique(database_id),
    groups_p = paste0(kidney_disease_group_short, ": ", hpo_id_group_p, collapse = " | "),
    .groups = "keep") %>%
  ungroup() %>%
  group_by(approved_symbol) %>%
  summarise(approved_symbol = unique(approved_symbol),
    groups_p = paste0(database_id, " (", groups_p, ")", collapse = "; "),
    .groups = "keep") %>%
  ungroup()

# TODO: check if it could be helpful to do this also per gene
# TODO: workflow: if the respective gene/entity from our kidney list is in the ClinGen curated list apply this group, if not apply the group with the highest score after reviewing the groups manually
############################################


############################################
# annotate pediatric vs adult onset (OMIM: HPO terms Adult onset HP:0003581, Pediatric onset HP:0410280, and children of terms)
# TODO: describe the logic of the analysis in comments step by step

# filter OMIM data for all the disease database ids with HPO terms for adult onset
phenotype_hpoa_filter_adult <- phenotype_hpoa %>%
   filter(hpo_id %in% hpo_list_adult$term) %>%
   select(database_id, hpo_id) %>%
   unique()

# annotate the OMIM data with the gene symbols
hpo_gene_list_adult <- phenotype_hpoa_filter_adult %>%
  left_join(omim_genemap2_disease_and_gene, by = c("database_id" = "disease_ontology_id"), relationship = "many-to-many") %>%
  filter(!is.na(approved_symbol)) %>%
  unique() %>%
  mutate(onset_group = "adult")

# filter OMIM data for all the disease database ids with HPO terms for non-adult onset
phenotype_hpoa_filter_non_adult <- phenotype_hpoa %>%
   filter(hpo_id %in% hpo_list_non_adult$term) %>%
   select(database_id, hpo_id) %>%
   unique()

# annotate the OMIM data with the gene symbols
hpo_gene_list_non_adult <- phenotype_hpoa_filter_non_adult %>%
  left_join(omim_genemap2_disease_and_gene, by = c("database_id" = "disease_ontology_id"), relationship = "many-to-many") %>%
  filter(!is.na(approved_symbol)) %>%
  unique() %>%
  mutate(onset_group = "non_adult")

# bind the two lists and compute the relative frequency of all onset related HPO terms in the OMIM entries associated with at least one of the terms per onset group
onset_hpo <- bind_rows(hpo_gene_list_adult, hpo_gene_list_non_adult) %>%
  select(onset_group, hpo_id) %>%
  group_by(onset_group) %>%
  mutate(hpo_id_count = n()) %>%
  ungroup() %>%
  group_by(onset_group, hpo_id) %>%
  summarise(hpo_id = unique(hpo_id),
    hpo_id_count = unique(hpo_id_count),
    hpo_id_group_count = n(),
    .groups = "keep") %>%
  ungroup() %>%
  mutate(hpo_id_group_frac = hpo_id_group_count / hpo_id_count) %>%
  select(-hpo_id_count, -hpo_id_group_count)

# now compute the sum of the relative frequency of the HPO terms in the two lists per gene and disease
# this will be the hpo id group p-value
# TODO: check if it could be helpful to do this also per gene
hpo_gene_list_onset_groups <- bind_rows(hpo_gene_list_adult, hpo_gene_list_non_adult) %>%
   select(-onset_group) %>%
  left_join(onset_hpo, by = c("hpo_id"), relationship = "many-to-many") %>%
  select(approved_symbol, database_id, onset_group, hpo_id_group_frac) %>%
  group_by(approved_symbol, database_id, onset_group) %>%
  summarise(approved_symbol = unique(approved_symbol),
    database_id = unique(database_id),
    onset_group = unique(onset_group),
    hpo_id_group_p = round(sum(hpo_id_group_frac), 3),
    .groups = "keep") %>%
  ungroup() %>%
  arrange(approved_symbol, desc(hpo_id_group_p))

# now summarize the results per gene
hpo_gene_list_onset_groups_summarized_for_join <- hpo_gene_list_onset_groups %>%
  group_by(approved_symbol, database_id) %>%
  summarise(approved_symbol = unique(approved_symbol),
    database_id = unique(database_id),
    groups_p = paste0(onset_group, ": ", hpo_id_group_p, collapse = " | "),
    .groups = "keep") %>%
  ungroup() %>%
  group_by(approved_symbol) %>%
  summarise(approved_symbol = unique(approved_symbol),
    groups_p = paste0(database_id, " (", groups_p, ")", collapse = "; "),
    .groups = "keep") %>%
  ungroup()
############################################


############################################
# syndromic vs non-syndromic (categories in OMIM: GROWTH, SKELETAL, NEUROLOGIC, HEAD & NECK; exclude: CARDIOVASCULAR, ABDOMEN, GENITOURINARY)
# TODO: there could be an error in the logic how to compute the relative frequency of the HPO terms in the two lists, this needs to be fixed
# TODO: Comment: I think it might still b e right, we just need to clearly explain what the value means and why the syndromic value differs from the sum of the other three values
phenotype_hpoa_filter_syndromic <- phenotype_hpoa %>%
   filter(hpo_id %in% hpo_list_syndromic$term) %>%
   select(database_id, hpo_id) %>%
   unique()

hpo_gene_list_syndromic <- phenotype_hpoa_filter_syndromic %>%
  left_join(omim_genemap2_disease_and_gene, by = c("database_id" = "disease_ontology_id"), relationship = "many-to-many") %>%
  filter(!is.na(approved_symbol)) %>%
  unique() %>%
  mutate(symptome_group = "syndromic")

# list to differentiate into the three organ systems
phenotype_hpoa_filter_growth <- phenotype_hpoa %>%
   filter(hpo_id %in% all_hpo_children_list_growth_term$term) %>%
   select(database_id, hpo_id) %>%
   unique()

hpo_gene_list_growth <- phenotype_hpoa_filter_growth %>%
  left_join(omim_genemap2_disease_and_gene, by = c("database_id" = "disease_ontology_id"), relationship = "many-to-many") %>%
  filter(!is.na(approved_symbol)) %>%
  unique() %>%
  mutate(symptome_group = "growth")

phenotype_hpoa_filter_skeletal <- phenotype_hpoa %>%
   filter(hpo_id %in% all_hpo_children_list_skeletal_term$term) %>%
   select(database_id, hpo_id) %>%
   unique()

hpo_gene_list_skeletal <- phenotype_hpoa_filter_skeletal %>%
  left_join(omim_genemap2_disease_and_gene, by = c("database_id" = "disease_ontology_id"), relationship = "many-to-many") %>%
  filter(!is.na(approved_symbol)) %>%
  unique() %>%
  mutate(symptome_group = "skeletal")

phenotype_hpoa_filter_neurologic <- phenotype_hpoa %>%
   filter(hpo_id %in% all_hpo_children_list_neurologic_term$term) %>%
   select(database_id, hpo_id) %>%
   unique()

hpo_gene_list_neurologic <- phenotype_hpoa_filter_neurologic %>%
  left_join(omim_genemap2_disease_and_gene, by = c("database_id" = "disease_ontology_id"), relationship = "many-to-many") %>%
  filter(!is.na(approved_symbol)) %>%
  unique() %>%
  mutate(symptome_group = "neurologic")

# bind the four lists and compute the relative frequency of the HPO terms in the three lists
symptome_hpo <- bind_rows(hpo_gene_list_syndromic, hpo_gene_list_growth, hpo_gene_list_skeletal, hpo_gene_list_neurologic) %>%
  select(symptome_group, hpo_id) %>%
  group_by(symptome_group) %>%
  mutate(hpo_id_count = n()) %>%
  ungroup() %>%
  group_by(symptome_group, hpo_id) %>%
  summarise(hpo_id = unique(hpo_id),
    hpo_id_count = unique(hpo_id_count),
    hpo_id_group_count = n(),
    .groups = "keep") %>%
  ungroup() %>%
  mutate(hpo_id_group_frac = hpo_id_group_count / hpo_id_count) %>%
  select(-hpo_id_count, -hpo_id_group_count)

hpo_gene_list_symptome_groups <- hpo_gene_list_syndromic %>%
   select(-symptome_group) %>%
  left_join(symptome_hpo, by = c("hpo_id"), relationship = "many-to-many") %>%
  select(approved_symbol, database_id, symptome_group, hpo_id_group_frac) %>%
  group_by(approved_symbol, database_id, symptome_group) %>%
  summarise(approved_symbol = unique(approved_symbol),
    database_id = unique(database_id),
    symptome_group = unique(symptome_group),
    hpo_id_group_p = round(sum(hpo_id_group_frac), 3),
    .groups = "keep") %>%
  ungroup() %>%
  arrange(approved_symbol, desc(hpo_id_group_p))

# summarize the results per gene
# always put syndromic first, then growth, then skeletal, then neurologic
hpo_gene_list_symptome_groups_summarized_for_join <- hpo_gene_list_symptome_groups %>%
  group_by(approved_symbol, database_id) %>%
  mutate(symptome_group = factor(symptome_group, levels = c("syndromic", "growth", "skeletal", "neurologic"))) %>%
  arrange(symptome_group) %>%
  summarise(approved_symbol = unique(approved_symbol),
    database_id = unique(database_id),
    groups_p = paste0(symptome_group, ": ", hpo_id_group_p, collapse = " | "),
    .groups = "keep") %>%
  ungroup() %>%
  group_by(approved_symbol) %>%
  summarise(approved_symbol = unique(approved_symbol),
    groups_p = paste0(database_id, " (", groups_p, ")", collapse = "; "),
    .groups = "keep") %>%
  ungroup()

# TODO: check if we need to do this analysis per entity instead per gene and then aggregate/ summarize per gene
############################################


############################################
# TODO: annotate MGI mouse phenotypes kidney
# use https://www.mousemine.org/mousemine/begin.do
# https://www.mousemine.org/mousemine/

# use R package: https://www.bioconductor.org/packages/release/bioc/vignettes/InterMineR/inst/doc/InterMineR.html
# TODO: need to find a table to show the genes names
# TODO: need to find a table to map the mouse phenotypes to HPO terms
############################################


############################################
# TODO: annotate StringDB interactions with strong evidence kidney genes
# use download tables: https://string-db.org/cgi/download?sessionId=bVJC28yKCRz5&species_text=Homo+sapiens
# eg physical: https://stringdb-downloads.org/download/protein.physical.links.v12.0/9606.protein.physical.links.v12.0.txt.gz
# sum the interactions for each gene with the high evidence list, the percentile normalize, additionally show a sting with interacting genes and scores
############################################


############################################
# TODO: move this to HGNC table
# TODO: annotate GTEx kidney expression
# use https://gtexportal.org/api/v2/redoc#tag/GTEx-Portal-API-Info
# see GitHub issue for cutoffs
# TODO: maybe add expression in embryonic kidney from somewhere (e.g. https://descartes.brotmanbaty.org/)
# descartes: https://atlas.fredhutch.org/data/bbi/descartes/human_gtex/tables/tissue_tpm/kidney.csv (this one for TPM)
# descartes: https://atlas.fredhutch.org/data/bbi/descartes/human_gtex/tables/tissue_percentage/kidney.csv
# descartes: https://atlas.fredhutch.org/data/bbi/descartes/human_gtex/counts/gene2ens.csv (for identifier mapping)

############################################


############################################
# TODO: add a column in the final table for publication (screening = 1 point, first clinical description = 2 points, clinical replication = 3 points)
# TODO: scoring logic: if only screening publication then the category cant be more then Limited, if there is a clinical description then the category can be Moderate, if there is a clinical replication then the category can be Definitive
# TODO: scoring logic: we further use the category "no known relation" for genes that are not trustworthy associated with kidney disease
# TODO: maybe annotate with publications (from OMIM) and GeneReviews

############################################


############################################
## save results
creation_date <- strftime(as.POSIXlt(Sys.time(),
  "UTC",
  "%Y-%m-%dT%H:%M:%S"), "%Y-%m-%d")

write_csv(results_genes_wider,
  file = paste0("merged/KidneyGenetics_AnnotateMergedTable.",
    creation_date,
    ".csv"),
  na = "NULL")

gzip(paste0("results/KidneyGenetics_AnnotateMergedTable.", creation_date, ".csv"),
  overwrite = TRUE)
############################################