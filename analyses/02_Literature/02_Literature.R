############################################
## load libraries
library(readxl) # needed to load Excel files
library(tidyverse) # needed for data processing
library(jsonlite) # needed for HGNC functions
library(officer) # needed for docx files
library(pdftools) # needed for pdf files
library(config) # needed for config loading
############################################


############################################
## define relative script path
project_topic <- "nephrology"
project_name <- "kidney-genetics"
script_path <- "/analyses/02_Literature/"

## read configs
config_vars <- config::get(file = Sys.getenv("CONFIG_FILE"),
    config = "default")
config_vars_path <- config::get(file = Sys.getenv("CONFIG_FILE"),
    config = project_topic)

## set working directory
setwd(paste0(config_vars_path$projectsdir, project_name, script_path))

## set global options
options(scipen = 999)
############################################


############################################
# load global functions
# hgnc functions
source("../functions/hgnc-functions.R", local = TRUE)
############################################


############################################
## load manually curated overview Excel table with download links and filter
pub_file <- "data/230220_Kidney_Genes_Publication_List.xlsx"

kidney_genes_publication_list <- read_excel(pub_file, skip = 4, na = "NA") %>%
  filter(!is.na(Type))
############################################


############################################
## download all files
# TODO: replace wininet method
# TODO: make download of supplements more stable (2,7,15,16), maybe get links from website
kidney_genes_publication_list_download <- kidney_genes_publication_list %>%
  rowwise() %>%
  mutate(downloaded = download.file(Download_link, paste0("data/downloads/PMID_", PMID, ".", Type), mode = "wb", quiet = TRUE, method = "wininet"))
############################################


############################################
## load and normalize the different files

###############
# Publication number 1: for PMID 35325889
number1_pmid_35325889_genes <- read_docx("data/PMID_35325889.docx") %>%
  docx_summary() %>%
  filter(content_type == "paragraph") %>%
  filter(str_detect(text, "ABCC6")) %>%
  select(gene = text) %>%
  separate_rows(gene, sep = ", ") %>%
  unique() %>%
  arrange(gene) %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(gene)) %>%
  mutate(source = "PMID_35325889") %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  select(approved_symbol, hgnc_id, gene_name_reported = gene, source)


###############
# Publication number 2: for PMID 34264297
number2_pmid_34264297_genes <- unzip("data/PMID_34264297.zip", list = TRUE, exdir = "data", overwrite = TRUE) %>%
  rowwise() %>%
  mutate(genes = list(read_excel(unzip("data/PMID_34264297.zip", files = Name, exdir = "data"), skip = 2) %>%
  select(Gene))) %>%
  ungroup() %>%
  unnest(genes) %>%
  select(gene = Gene) %>%
  mutate(gene = str_replace_all(gene, "orf", "ORF")) %>%
  mutate(gene = str_remove_all(gene, "[ \\/].+")) %>%
  mutate(gene = str_remove_all(gene, "[\\(\\*].+")) %>%
  mutate(gene = str_remove_all(gene, "\\*")) %>%
  filter(gene != "") %>%
  filter(!is.na(gene)) %>%
  filter(!str_detect(gene, "[:-]")) %>%
  filter(!str_detect(gene, "[a-z,]")) %>%
  unique() %>%
  arrange(gene) %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(gene)) %>%
  mutate(source = "PMID_34264297") %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  select(approved_symbol, hgnc_id, gene_name_reported = gene, source)


###############
# Publication number 3: for PMID 36035137
number3_pmid_36035137_genes <- read_excel("data/PMID_36035137.xlsx", skip = 1) %>%
  select(gene = Gene) %>%
  unique() %>%
  arrange(gene) %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(gene)) %>%
  mutate(source = "PMID_36035137") %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  select(approved_symbol, hgnc_id, gene_name_reported = gene, source)


###############
# Publication number 5:
number5_pmid_33664247_genes <- pdf_text("data/PMID_33664247.pdf") %>%
  read_lines(skip = 15) %>%
  str_squish() %>%
  as_tibble() %>%
  select(gene = value) %>%
  separate_rows(gene, sep = ", | |,") %>%
  filter(gene != "") %>%
  filter(!(gene %in% c("1", "1.", "2", "3", "552"))) %>%
  filter(!str_detect(gene, pattern = "Panel|Gene|ADTKD|aHUS|C3|GN|Alport|syndrome|ARPKD|BORS|CAKUT|Cystinosis|Nephronophthisis|\\&|related|ciliopathies|Nephrotic|Tubulopathies|list|23|PLG|removed|EVC|Supplementary|Figure|Distribution|mutation|types|of|variants|uncertain|significance|VOUS|probands|Variant|classification|was|based|2015|ACMG|guidelines|Abbreviated|mutation|types|are|as|follows|CNV|copy|number|variation|indels|insertions|deletions|or|\\(n|\\=")) %>%
  mutate(gene = str_remove(gene, "_.+")) %>%
  unique() %>%
  arrange(gene) %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(gene)) %>%
  mutate(source = "PMID_33664247") %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  select(approved_symbol, hgnc_id, gene_name_reported = gene, source)


###############
# Publication number 7: for PMID 30476936
number7_pmid_30476936_genes <- read_excel("data/PMID_30476936.xlsx") %>%
  select(gene = `Gene Name`) %>%
  unique() %>%
  arrange(gene) %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(gene)) %>%
  mutate(source = "PMID_30476936") %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  select(approved_symbol, hgnc_id, gene_name_reported = gene, source)


###############
# Publication number 9: for PMID 31509055
number9_pmid_31509055_genes <- pdf_text("data/PMID_31509055.pdf") %>%
  read_lines() %>%
  str_squish() %>%
  as_tibble() %>%
  select(gene = value) %>%
  separate_rows(gene, sep = ", | |,") %>%
  filter(gene != "") %>%
  filter(!str_detect(gene, pattern = "Supplementary|Table|1:|List|of|genes|included|in|gene|panel|\\(n|\\=")) %>%
  mutate(gene = str_remove(gene, "_.+")) %>%
  unique() %>%
  arrange(gene) %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(gene)) %>%
  mutate(source = "PMID_31509055") %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  select(approved_symbol, hgnc_id, gene_name_reported = gene, source)


###############
# Publication number 10: for PMID 31822006
number10_pmid_31822006_genes <- pdf_text("data/PMID_31822006.pdf")[1:9] %>%
  read_lines(skip = 31) %>%
  str_squish() %>%
  as_tibble() %>%
  select(gene = value) %>%
  mutate(gene = str_remove_all(gene, " .+")) %>%
  filter(gene != "") %>%
  filter(!str_detect(gene, "[a-z,]")) %>%
  filter(!str_detect(gene, "^[0-9]$")) %>%
  unique() %>%
  arrange(gene) %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(gene)) %>%
  mutate(source = "PMID_31822006") %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  select(approved_symbol, hgnc_id, gene_name_reported = gene, source)


###############
# Publication number 12: for PMID 29801666
number12_pmid_29801666_genes <- read_docx("data/PMID_29801666.docx") %>%
  docx_summary() %>%
  filter(cell_id == 1) %>%
  filter(text != "Gene") %>%
  select(gene = text) %>%
  unique() %>%
  arrange(gene) %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(gene)) %>%
  mutate(source = "PMID_29801666") %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  select(approved_symbol, hgnc_id, gene_name_reported = gene, source)


###############
# Publication number 13: for PMID 31027891
number13_pmid_31027891_genes <- read_docx("data/PMID_31027891.docx") %>%
  docx_summary() %>%
  filter(cell_id == 8) %>%
  filter(text != "GENE SYMBOL") %>%
  select(gene = text) %>%
  separate_rows(gene, sep = ", ") %>%
  unique() %>%
  arrange(gene) %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(gene)) %>%
  mutate(source = "PMID_31027891") %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  select(approved_symbol, hgnc_id, gene_name_reported = gene, source)


###############
# Publication number 14: for PMID 26862157
number14_pmid_26862157_genes <- pdf_text("data/PMID_26862157.pdf")[1:4] %>%
  read_lines(skip = 9) %>%
  str_squish() %>%
  as_tibble() %>%
  select(gene = value) %>%
  mutate(gene = str_remove_all(gene, " .+")) %>%
  mutate(gene = str_replace_all(gene, "orf", "ORF")) %>%
  filter(gene != "") %>%
  filter(!str_detect(gene, "[a-z,]")) %>%
  filter(!str_detect(gene, "^[0-9]+$")) %>%
  unique() %>%
  arrange(gene) %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(gene)) %>%
  mutate(source = "PMID_26862157") %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  select(approved_symbol, hgnc_id, gene_name_reported = gene, source)


###############
# Publication number 15: for PMID 33532864
number15_pmid_33532864_genes <- read_docx("data/PMID_33532864.docx") %>%
  docx_summary() %>%
  filter(cell_id == 1) %>%
  filter(text != "Gene") %>%
  select(gene = text) %>%
  separate_rows(gene, sep = ", ") %>%
  slice_head(n = 316) %>%
  unique() %>%
  arrange(gene) %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(gene)) %>%
  mutate(source = "PMID_33532864") %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  select(approved_symbol, hgnc_id, gene_name_reported = gene, source)


###############
# Publication number 16: for PMID 35005812
number16_pmid_35005812_genes <- read_excel("data/PMID_35005812.xlsx", skip = 1, na = "NA") %>%
  select(gene = `...2`) %>%
  mutate(gene = str_remove_all(gene, "[\\( ].+")) %>%
  filter(!is.na(gene)) %>%
  unique() %>%
  arrange(gene) %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(gene)) %>%
  mutate(source = "PMID_35005812") %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  select(approved_symbol, hgnc_id, gene_name_reported = gene, source)

############################################


############################################
# Load all publication results into one tibble
literature_genes <- kidney_genes_publication_list %>%
  select(Number, PMID) %>%
  rowwise() %>%
  mutate(result = list(eval(parse(text = paste0("number", Number, "_pmid_", PMID, "_genes"))))) %>%
  unnest(result) %>%
  select(approved_symbol, hgnc_id, gene_name_reported, source) %>%
  group_by(approved_symbol, hgnc_id) %>%
  summarise(gene_name_reported = paste(unique(gene_name_reported), collapse = " | "),
    source = paste(unique(source), collapse = " | "),
    publication_count = n(),
    .groups = "keep") %>%
  ungroup() %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(approved_symbol)) %>%
  mutate(at_least_two_publications = (publication_count > 1)) %>%
  select(approved_symbol,
    hgnc_id,
    gene_name_reported,
    source,
    source_count = publication_count,
    source_evidence = at_least_two_publications)

# TODO: normalize PMID representation

# TODO: normalize source_evidence to 0/1 as percentiles
# TODO: write a function for this normalization step

############################################


############################################
## save results
# TODO: gzip csv result files
# TODO: add kidney_genes_publication_list_download with checksums to results
creation_date <- strftime(as.POSIXlt(Sys.time(),
  "UTC",
  "%Y-%m-%dT%H:%M:%S"), "%Y-%m-%d")

write_csv(literature_genes,
  file = paste0("results/02_Literature_genes.",
    creation_date,
    ".csv"),
  na = "NULL")
############################################