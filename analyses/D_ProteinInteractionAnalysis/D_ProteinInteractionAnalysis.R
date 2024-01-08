############################################
## load libraries
library(tidyverse)
library(STRINGdb)
library(R.utils)
library(plotly)
library(network)
library(ggplot2)
library(GGally)
############################################

# TODO: set in config?
string_db_version <- "12.0"
string_db_files_path <- "../shared/"
min_gene_number_per_cluster <- 60  # TODO: StackError if min_gene_number_per_cluster < 60

############################################
## define relative script path
project_topic <- "nephrology"
project_name <- "kidney-genetics"
script_path <- "/analyses/D_ProteinInteractionAnalysis/"

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
source("../functions/file-functions.R", local = TRUE)
# helper functions
source("../functions/helper-functions.R", local = TRUE)
# cluster analysis functions
source("../functions/protein_interaction_analysis-functions.R", local = TRUE)
############################################


############################################
# compute date only once or somehow in config
current_date <- get_current_date_iso8601()
############################################


############################################
# create a directory to save the plots of the analysis
output_dir <- "results"
dir.create(output_dir, showWarnings = TRUE)
############################################


############################################
# load annotated HGNC table, select required columns
hgnc_annotated_path <- get_newest_file(file_basename = "non_alt_loci_set_coordinates",
                                       folder = "../B_AnnotationHGNC/results")

hgnc_annotated <- read.csv(gzfile(hgnc_annotated_path)) %>% 
  dplyr::select(hgnc_id, symbol, STRING_id) %>% 
  mutate(hgnc_id = as.numeric(sub("HGNC:", "", hgnc_id)))

# load annotated high evidence genes, spread hpo_id_group_p 
high_evidence_annotated_path <- get_newest_file(file_basename = "C_high_evidence_annotated_csv_table",
                                                folder = "../C_AnnotateMergedTable/results")

high_evidence_annotated <- read.csv(gzfile(high_evidence_annotated_path), na.strings = c("NULL", NaN)) %>% 
  dplyr::select(hgnc_id, approved_symbol, clinical_groups_p) %>% 
  mutate(clinical_groups_p = sub(".*\\((.*)\\).*", "\\1", clinical_groups_p)) %>% 
  separate_rows(clinical_groups_p, sep = "\\|") %>% 
  filter(!is.na(clinical_groups_p)) %>% 
  separate(clinical_groups_p, into = c("kidney_disease_group_short", "hpo_id_group_p"), sep = ":\\s", remove = TRUE) %>% 
  mutate(kidney_disease_group_short = str_replace_all(kidney_disease_group_short, "\\s+", ""),
         hpo_id_group_p = as.numeric(hpo_id_group_p))

# join annotated high evidence genes with STRING id
kid_groups <- high_evidence_annotated %>% 
  left_join(hgnc_annotated, by = c("hgnc_id")) %>% 
  filter(!is.na(STRING_id))
############################################


############################################
## Cluster analysis with FULL network
# download STRING db protein links 
protein_links_url <- paste0("https://stringdb-downloads.org/download/protein.links.v", string_db_version, "/9606.protein.links.v", string_db_version, ".txt.gz")
protein_links_file <- paste0("../shared/9606.protein.links.v", string_db_version, ".txt.gz")

if (!file.exists(protein_links_file)) {
  download.file(protein_links_url, destfile = protein_links_file, method = "auto")
  message("STRING protein links downloaded successfully.")
} else {
  message(paste0("\"", protein_links_file, "\" already exists."))
}

# instantiate new STRING db reference class 
string_db_full <- STRINGdb::STRINGdb$new(version = string_db_version,
                                    species = 9606,
                                    score_threshold = 100,
                                    input_directory = string_db_files_path)

# get STRING IDs of genes that have a kidney disease group
STRING_id_vec <- unique(kid_groups$STRING_id)

# get list of sublists that contain STRING subclusters 
STRING_clusters_list <- get_STRING_clusters(STRING_id_vec = STRING_id_vec, min_number = min_gene_number_per_cluster)

# create a df that contains full index of each gene with the subclusters
cluster_index_df <- STRING_clusters_list %>% 
  unlist() %>% 
  as.data.frame() %>%
  setNames("STRING_id") %>% rowwise() %>% 
  mutate(cluster_index = get_full_index(STRING_clusters_list, STRING_id))

# get the maximum cluster depth
max_depth <- max(sapply(strsplit(cluster_index_df$cluster_index, "-"), length))

# separate 'cluster_index' into multiple columns
cluster_index_df <- cluster_index_df %>% 
  mutate(split_col = cluster_index) %>% 
  separate(split_col, into = paste0("hierarchy_", 1:max_depth), sep = "-", fill = 'right')

# join with gene symbol and HGNC ID
cluster_index_df <- cluster_index_df %>% 
  left_join(hgnc_annotated[c("STRING_id", "symbol", "hgnc_id")], by = "STRING_id")

# write csv
write_csv(cluster_index_df,
          file = paste0("results/STRING_cluster_indices_min_gene_number-", min_gene_number_per_cluster, "-", current_date, ".csv"))

gzip(paste0("results/STRING_cluster_indices_min_gene_number-", min_gene_number_per_cluster, "-", current_date, ".csv"),
     overwrite = TRUE)
############################################


############################################
# Plot distribution of kidney disease groups within subcluster

# get most probable kidney disease group per gene and join with cluster index df
max_groups_full <- kid_groups %>%
  group_by(approved_symbol) %>%
  filter(hpo_id_group_p == max(hpo_id_group_p)) %>% 
  left_join(cluster_index_df[c("hgnc_id", "cluster_index")], by = "hgnc_id") %>% 
  filter(!is.na(cluster_index))
# NOTE: if genes have more than one group with same and highest hpo_id_group_p => keep all groups 

# write csv
write_csv(max_groups_full,
          file = paste0("results/disease_group_STRING_cluster_indices_min_gene_number-", min_gene_number_per_cluster, "-", current_date, ".csv"))

gzip(paste0("results/disease_group_STRING_cluster_indices_min_gene_number-", min_gene_number_per_cluster, "-", current_date, ".csv"),
     overwrite = TRUE)

# example plot of kidney disease group distribution within subcluster
ex1 <- plot_disease_group_distribution(subcluster="3-1", max_groups_full)
############################################


############################################
# Plot interaction network


plot_network_of_index_gene <- function(index_gene, string_db, min_comb_score, STRING_id_vec, disease_group_df){
  
  # get all STRING interactions between all genes in STRING_id_vec with a minimum combined score in STRING
  all_interactions <- get_all_interactions_above_score(STRING_id_vec = STRING_id_vec,
                                                       string_db = string_db,
                                                       min_comb_score = min_comb_score)
  # get all contacts of the index_gene above the minimum combined score in STRING
  
  all_contacts_to_index <- get_all_contacts(index_gene = index_gene,
                                            connection_df = all_interactions,
                                            min_comb_score = min_comb_score)
  
  # create an edgelist of all contacts of the index gene (direct/indirect contacts)
  edgelist <- create_edgelist(connection_df = all_interactions, 
                              all_contacts_to_index = all_contacts_to_index, 
                              min_comb_score = min_comb_score, 
                              symbol_annotation_df = distinct(disease_group_df[c("STRING_id", "symbol")]))
  
  
  # plot interaction network of index gene
  interaction_plot <- plot_interaction_network(edgelist = edgelist,
                                               disease_group_df = disease_group_df)
  
  return(interaction_plot)
}

# 
# # example network plot for "MAPK1" ("9606.ENSP00000215832")
# ia_plot <- plot_network_of_index_gene(index_gene = "9606.ENSP00000215832", 
#                            string_db = string_db_full, 
#                            min_comb_score = 980, 
#                            STRING_id_vec = STRING_id_vec, 
#                            disease_group_df = max_groups_full)
# 
# ia_plot
