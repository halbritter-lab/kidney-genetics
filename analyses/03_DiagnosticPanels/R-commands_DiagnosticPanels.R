############################################
## load libraries
library(readr)
library(readxl)
library(tidyverse)
library(rvest)
library(jsonlite)
library(curl)
library(httr)
library(config) # needed for config loading
library(webdriver) # needed for headless browsing
############################################


############################################
## define relativescript path
script_path <- "kidney-genetics/analyses/03_DiagnosticPanels/"
## read config
config_vars <- config::get(file = Sys.getenv("CONFIG_FILE"))
## set working directory
setwd(paste0(config_vars$projectsdir, script_path))
## set global options
options(scipen = 999)
############################################


############################################
# load global functions
# hgnc functions
source("../functions/hgnc-functions.R", local = TRUE)
source("../functions/phantomjs-functions.R", local = TRUE)
source("../functions/natera-functions.R", local = TRUE)
source("../functions/blueprintgenetics-functions.R", local = TRUE)
############################################


############################################
## download web urls

# load the list of sources
diagnostic_panels_list <- read_excel("data/diagnostic_panels_list.xlsx") %>%
	filter(use == "yes")

# download using phantomJS
diagnostic_panels <- diagnostic_panels_list %>%
	rowwise() %>%
	mutate(filename_download = download_url_by_phantomjs(diagnostic_panel_source, diagnostic_panel_name, type)) %>%
	ungroup()
############################################


############################################
## compute standardized gene list per panel

## 1) centogene_nephrology
url <- (diagnostic_panels %>%
			filter(diagnostic_panel_name == "centogene_nephrology"))$filename_download

centogene_nephrology <- read_html(url)

centogene_nephrology_genes_nodes <- centogene_nephrology %>%
	html_nodes(xpath = '//*[@id="t3m-Modal--74334"]//table[1]') %>%
	html_table()
	
centogene_nephrology_genes <- centogene_nephrology_genes_nodes[[1]] %>%
	select(Genes = Gene) %>%
	unique() %>%
	mutate(Panel = "centogene_nephrology") %>%
	mutate(Source = url) %>%
	mutate(Genes = str_remove_all(Genes, "[\\*\\#]")) %>%
	mutate(hgnc_id = hgnc_id_from_symbol_grouped(Genes)) %>%
	mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
	select(approved_symbol, hgnc_id, gene_name_reported = Genes, source = Source, panel = Panel)


## 2) cegat_kidney_diseases
url <- (diagnostic_panels %>%
			filter(diagnostic_panel_name == "cegat_kidney_diseases"))$filename_download

cegat_kidney_diseases <- read_html(url)

cegat_kidney_diseases_genes <- cegat_kidney_diseases %>%
    html_nodes(xpath = '//h2[contains(text(),"Gene Directory")]//following::em') %>%
	html_text() %>%
	tibble(`Genes` = .) %>%
	separate_rows(., Genes, convert = TRUE) %>%
	mutate(Panel = "cegat_kidney_diseases_genes") %>%
	mutate(Source = url) %>%
	mutate(Genes = str_remove_all(Genes, "[\\*\\#]")) %>%
	mutate(hgnc_id = hgnc_id_from_symbol_grouped(Genes)) %>%
	mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
	select(approved_symbol, hgnc_id, gene_name_reported = Genes, source = Source, panel = Panel)


## 3) preventiongenetics_comprehensive_inherited_kidney_diseases_panel
url <- (diagnostic_panels %>%
			filter(diagnostic_panel_name == "preventiongenetics_comprehensive_inherited_kidney_diseases_panel"))$filename_download

preventiongenetics_comprehensive_inherited_kidney_diseases_panel <- read_html(url)

preventiongenetics_comprehensive_inherited_kidney_diseases_panel_genes <- preventiongenetics_comprehensive_inherited_kidney_diseases_panel %>%
	html_nodes(xpath = '//*[@id="genes-div"]//table[1]') %>%
	html_table()
	
preventiongenetics_comprehensive_inherited_kidney_diseases_panel_genes <- preventiongenetics_comprehensive_inherited_kidney_diseases_panel_genes[[1]] %>%
	select(Genes = `Official Gene Symbol`) %>%
	unique() %>%
	mutate(Panel = "preventiongenetics_comprehensive_inherited_kidney_diseases_panel") %>%
	mutate(Source = url) %>%
	mutate(Genes = str_remove_all(Genes, "[\\*\\#]")) %>%
	mutate(hgnc_id = hgnc_id_from_symbol_grouped(Genes)) %>%
	mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
	select(approved_symbol, hgnc_id, gene_name_reported = Genes, source = Source, panel = Panel)


## 4) invitae_progressive_renal_disease_panel
url <- (diagnostic_panels %>%
			filter(diagnostic_panel_name == "invitae_progressive_renal_disease_panel"))$filename_download

invitae_progressive_renal_disease_panel <- read_html(url)

invitae_progressive_renal_disease_panel_genes <- invitae_progressive_renal_disease_panel %>%
    html_nodes(xpath = '//meta[contains(@content,"ACE")]') %>% 
    html_attr("content") %>%
	tibble(`Genes` = .) %>%
	separate_rows(., Genes, convert = TRUE) %>%
	mutate(Panel = "invitae_progressive_renal_disease_panel") %>%
	mutate(Source = url) %>%
	mutate(Genes = str_remove_all(Genes, "[\\*\\#]")) %>%
	mutate(hgnc_id = hgnc_id_from_symbol_grouped(Genes)) %>%
	mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
	select(approved_symbol, hgnc_id, gene_name_reported = Genes, source = Source, panel = Panel)


## 5) invitae_expanded_renal_disease_panel
url <- (diagnostic_panels %>%
			filter(diagnostic_panel_name == "invitae_expanded_renal_disease_panel"))$filename_download

invitae_expanded_renal_disease_panel <- read_html(url)

invitae_expanded_renal_disease_panel_genes <- invitae_expanded_renal_disease_panel %>%
    html_nodes(xpath = '//meta[contains(@content,"ACE")]') %>% 
    html_attr("content") %>%
	tibble(`Genes` = .) %>%
	separate_rows(., Genes, convert = TRUE) %>%
	mutate(Panel = "invitae_expanded_renal_disease_panel") %>%
	mutate(Source = url) %>%
	mutate(Genes = str_remove_all(Genes, "[\\*\\#]")) %>%
	mutate(hgnc_id = hgnc_id_from_symbol_grouped(Genes)) %>%
	mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
	select(approved_symbol, hgnc_id, gene_name_reported = Genes, source = Source, panel = Panel)

## 6) mgz_nephrologie
url <- (diagnostic_panels %>%
			filter(diagnostic_panel_name == "mgz_nephrologie"))$filename_download

mgz_nephrologie <- read_html(url)

mgz_nephrologie_genes <- mgz_nephrologie %>%
    html_nodes(xpath = '//div[contains(@class,"panel_gene")]') %>%
	html_text() %>%
	str_remove_all(., "[\\n\\r ]+") %>%
	tibble(`Genes` = .) %>%
	mutate(Panel = "mgz_nephrologie") %>%
	mutate(Source = url) %>%
	mutate(Genes = str_remove_all(Genes, "[\\*\\#]")) %>%
	filter(!str_detect(Genes, "\\(")) %>%
	mutate(hgnc_id = hgnc_id_from_symbol_grouped(Genes)) %>%
	mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
	select(approved_symbol, hgnc_id, gene_name_reported = Genes, source = Source, panel = Panel) %>%
	unique()


## 7) mvz_nierenerkrankungen
# TODO: implement loading directly from downloaded file
url <- (diagnostic_panels %>%
			filter(diagnostic_panel_name == "mvz_nierenerkrankungen"))$diagnostic_panel_source

mvz_nierenerkrankungen <- fromJSON(url)

mvz_nierenerkrankungen_genes <- tibble(mvz_nierenerkrankungen$Gene) %>%
	select(`Genes` = Genname) %>%
	mutate(Panel = "mvz_nierenerkrankungen") %>%
	mutate(Source = url) %>%
	mutate(hgnc_id = hgnc_id_from_symbol_grouped(Genes)) %>%
	mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
	select(approved_symbol, hgnc_id, gene_name_reported = Genes, source = Source, panel = Panel) %>%
	unique()


## 8) natera_renasight_comprehensive_kidney_gene_panel
url <- (diagnostic_panels %>%
			filter(diagnostic_panel_name == "natera_renasight_comprehensive_kidney_gene_panel",
			method == "read_html"))$diagnostic_panel_source

api_url <- (diagnostic_panels %>%
			filter(diagnostic_panel_name == "natera_renasight_comprehensive_kidney_gene_panel",
			method == "curl"))$diagnostic_panel_source

natera_renasight_comprehensive_kidney_gene_panel <- list()
for(i in 1:natera_renasight_get_last_page_number()){
  natera_page <- natera_renasight_get_genes_from_page(i, api_url)
  natera_renasight_comprehensive_kidney_gene_panel[[i+1]] <- natera_page
}

natera_renasight_comprehensive_kidney_gene_panel_genes <- natera_renasight_comprehensive_kidney_gene_panel %>%
	unlist() %>%
	tibble(`Genes` = .) %>%
	mutate(Panel = "natera_renasight_comprehensive_kidney_gene_panel") %>%
	mutate(Source = url) %>%
	mutate(hgnc_id = hgnc_id_from_symbol_grouped(Genes)) %>%
	mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
	select(approved_symbol, hgnc_id, gene_name_reported = Genes, source = Source, panel = Panel)


## 9) mayocliniclabs_renal_genetics
url <- (diagnostic_panels %>%
			filter(diagnostic_panel_name == "mayocliniclabs_renal_genetics"))$filename_download
mayocliniclabs_renal_genetics <- read_html(url)

mayocliniclabs_renal_genetics_genes <- mayocliniclabs_renal_genetics %>%
    html_nodes(xpath = '//span[contains(text(),"This test utilizes next-generation sequencing to detect")]//i') %>%
	html_text() %>%
	str_remove_all(., "[\\n ]+") %>%
	tibble(`Genes` = .) %>%
	separate_rows(., Genes, convert = TRUE) %>%
	mutate(Panel = "mayocliniclabs_renal_genetics") %>%
	mutate(Source = url) %>%
	mutate(Genes = str_remove_all(Genes, "[\\*\\#\\.]")) %>%
	filter(!str_detect(Genes, "\\(")) %>%
	unique() %>%
	filter(Genes != "") %>%
	filter(Genes != "N27") %>%
	mutate(hgnc_id = hgnc_id_from_symbol_grouped(Genes)) %>%
	mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
	select(approved_symbol, hgnc_id, gene_name_reported = Genes, source = Source, panel = Panel)


## 10) blueprintgenetics_nephrology
# TODO: Does not work anymore, fix
# TODO: implement download of all subpanels as files
url <- (diagnostic_panels %>%
			filter(diagnostic_panel_name == "blueprintgenetics_nephrology"))$diagnostic_panel_source

# generate list of blueprintgenetics sub panels
blueprintgenetics_panel_list <- diagnostic_panels %>%
			filter(diagnostic_panel_name == "blueprintgenetics_nephrology") %>%
			select(diagnostic_panel_name, subpanel_name, subpanel_source) %>%
			separate_rows(c(subpanel_name, subpanel_source) , sep = ", ", convert = TRUE) %>%
			mutate(subpanel_name = paste0(diagnostic_panel_name, "_", subpanel_name)) %>%
			mutate(type = "html")

# download using phantomJS
blueprintgenetics_panel_panels <- blueprintgenetics_panel_list %>%
	rowwise() %>%
	mutate(filename_download = download_url_by_phantomjs(subpanel_source, subpanel_name, type)) %>%
	ungroup()

# read all blueprintgenetics nephrology sub panels
blueprintgenetics_panel_panels_genes <- blueprintgenetics_panel_panels %>%
	rowwise() %>%
	mutate(gene_list = list(blueprintgenetics_panel_extract_genes(filename_download))) %>%
	ungroup()

blueprintgenetics_nephrology_genes <- blueprintgenetics_panel_panels_genes$gene_list %>%
	unlist() %>%
	tibble(`Genes` = .) %>%
	mutate(Panel = "blueprintgenetics_nephrology") %>%
	mutate(Source = url) %>%
	mutate(Genes = str_remove_all(Genes, "[\\*\\#\\.]")) %>%
	unique() %>%
	mutate(hgnc_id = hgnc_id_from_symbol_grouped(Genes)) %>%
	mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
	select(approved_symbol, hgnc_id, gene_name_reported = Genes, source = Source, panel = Panel)
############################################


############################################
## bind all tables and summaize
## compute count of panels a gene is reported in
## compute diagostic panel evidence as genes reported in at least 2 panel sources
all_diagnostic_panels_genes <- bind_rows(centogene_nephrology_genes,
		cegat_kidney_diseases_genes,
		preventiongenetics_comprehensive_inherited_kidney_diseases_panel_genes,
		invitae_progressive_renal_disease_panel_genes,
		invitae_expanded_renal_disease_panel_genes,
		mgz_nephrologie_genes,
		mvz_nierenerkrankungen_genes,
		natera_renasight_comprehensive_kidney_gene_panel_genes,
		mayocliniclabs_renal_genetics_genes,
		blueprintgenetics_nephrology_genes) %>%
	group_by(approved_symbol) %>%
	summarise(panel_diagnostic = paste(unique(panel), collapse = "; "),
		panel_diagnostic_count = n(),
		gene_name_reported = paste(unique(gene_name_reported), collapse = " | ")) %>%
	ungroup() %>%
	mutate(hgnc_id = paste0("HGNC:", hgnc_id_from_symbol_grouped(approved_symbol))) %>%
	mutate(at_least_two_panels = (panel_diagnostic_count > 1))

all_diagnostic_panels_genes_format <- all_diagnostic_panels_genes %>%
  select(approved_symbol, hgnc_id, gene_name_reported, source = panel_diagnostic, source_count = panel_diagnostic_count, source_evidence = at_least_two_panels)
############################################


############################################
## save results
creation_date <- strftime(as.POSIXlt(Sys.time(), "UTC", "%Y-%m-%dT%H:%M:%S"), "%Y-%m-%d")
write_csv(all_panels_genes, file = paste0("results/03_DiagnosticPanels_genes.", creation_date, ".csv"), na="NULL")
write_csv(diagnostic_panels, file = paste0("results/03_DiagnosticPanels_list.", creation_date, ".csv"), na="NULL")
############################################
