############################################
## load libraries
library(tidyverse)
library(STRINGdb)
library(ggplot2)
############################################

# TODO: set in config?
string_db_version <- "12.0"
string_db_files_path <- "."

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
############################################


############################################
# compute date only once or somehow in config
current_date <- get_current_date_iso8601()
############################################


############################################
# create a directory to save the plots of the analysis
output_dir <- "plots"
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
## analysis with FULL network
# download STRING db protein links 
protein_links_url <- paste0("https://stringdb-downloads.org/download/protein.links.v", string_db_version, "/9606.protein.links.v", string_db_version, ".txt.gz")
protein_links_file <- paste0("9606.protein.links.v", string_db_version, ".txt.gz")

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

# get clusters
clusters_list_full <- string_db_full$get_clusters(unique(kid_groups$STRING_id),
                                        algorithm = "walktrap") # NOTE: some STRING_ids get lost here => no group for them???

clusters_tibble_full <- tibble(clusters_list_full) %>% 
  select(STRING_id = clusters_list_full) %>% 
  mutate(cluster = row_number()) %>% 
  unnest_longer(col = "STRING_id")

max_groups_full <- kid_groups %>%
  group_by(approved_symbol) %>%
  filter(hpo_id_group_p == max(hpo_id_group_p)) %>% 
  left_join(clusters_tibble_full, by = "STRING_id") %>% 
  filter(!is.na(cluster))
# NOTE: if genes have more than one group with same and highest hpo_id_group_p => keep all groups 

# group the data by 'cluster' and calculate the counts for each group
cluster_counts <- table(max_groups_full$cluster)

# create an empty list to store the pie charts
pie_charts <- list()

# loop through each cluster group and create a pie chart
for (i in unique(max_groups_full$cluster)) {
  # filter the data for the current cluster group
  subset_df <- subset(max_groups_full, cluster == i)
  
  # create a pie chart for the current cluster group
  pie_chart <- ggplot(subset_df, aes(x = "", fill = kidney_disease_group_short)) +
    geom_bar(width = 1) +
    coord_polar(theta = "y") +
    labs(title = paste("Cluster", i))
  
  # add the pie chart to the list
  pie_charts[[i]] <- pie_chart
}

# calculate the relative frequency of each cluster
cluster_freq <- prop.table(table(max_groups_full$cluster))

# define a custom color palette for kidney_disease_group_short
custom_colors <- c("tubulopathy" = "Red", 
                   "glomerulopathy" = "Green", 
                   "cancer" = "Blue", 
                   "cakut" = "Purple", 
                   "cyst_cilio" = "Orange", 
                   "complement" = "Yellow",
                   "nephrocalcinosis" = "Grey")  # Add more colors if needed

# create a data frame for pie chart creation
pie_data <- data.frame(cluster = names(cluster_freq), freq = cluster_freq)

# create and save a pie chart for each cluster group
for (i in seq_along(pie_data$cluster)) {
  subset_df <- subset(max_groups_full, cluster == pie_data$cluster[i])
  
  # calculate the total number of instances in the current cluster
  total_instances <- nrow(subset_df)
  
  pie_chart <- ggplot(subset_df, aes(x = "", fill = kidney_disease_group_short)) +
    geom_bar(width = 1) +
    scale_fill_manual(values = custom_colors) +
    coord_polar(theta = "y") +
    labs(title = paste("Full network - Cluster", pie_data$cluster[i], "- Total Instances:", total_instances))
  
  # define the filename based on cluster number
  filename <- paste0(output_dir, "/protein_interactions_full_network_cluster_", pie_data$cluster[i], "_pie_chart.", current_date, ".png")
  
  # save the pie chart as a PNG file
  ggsave(filename, plot = pie_chart, width = 5, height = 5)  
}

# TODO: get enrichment pathways etc
# group1 <- max_groups_full %>% filter(cluster == 1)
# 
# string_db_full$get_enrichment(group1$STRING_id) %>%
#   tibble() %>% View
# junk <- string_db_full$get_enrichment(group1$STRING_id) %>%
#   tibble()


############################################
## Analysis with PHYSICAL network
# download STRING db protein physical links 
protein_phys_links_url <- paste0("https://stringdb-static.org/download/protein.physical.links.v", string_db_version, "/9606.protein.physical.links.v", string_db_version, ".txt.gz")
protein_phys_links_file <- paste0("9606.protein.physical.links.v", string_db_version, ".txt.gz")

if (!file.exists(protein_phys_links_file)) {
  download.file(protein_phys_links_url, destfile = protein_phys_links_file, method = "auto")
  message("STRING protein physical links downloaded successfully.")
} else {
  message(paste0("\"", protein_phys_links_file, "\" already exists."))
}

# instantiate new STRING db reference class 
string_db_phys <- STRINGdb::STRINGdb$new(version = string_db_version,
                                         species = 9606,
                                         score_threshold = 100,
                                         network_type = "physical",
                                         input_directory = string_db_files_path)

# get clusters
clusters_list_phys <- string_db_phys$get_clusters(unique(kid_groups$STRING_id),
                                                  algorithm = "walktrap") # NOTE: some STRING_ids get lost here => no group for them???

clusters_tibble_phys <- tibble(clusters_list_phys) %>% 
  select(STRING_id = clusters_list_phys) %>% 
  mutate(cluster = row_number()) %>% 
  unnest_longer(col = "STRING_id")

max_groups_phys <- kid_groups %>%
  group_by(approved_symbol) %>%
  filter(hpo_id_group_p == max(hpo_id_group_p)) %>% 
  left_join(clusters_tibble_phys, by = "STRING_id") %>% 
  filter(!is.na(cluster))
# NOTE: if genes have more than one group with same and highest hpo_id_group_p => keep all groups 

# group the data by 'cluster' and calculate the counts for each group
cluster_counts <- max_groups_phys %>% 
  group_by(cluster) %>% 
  summarise(count = n())

# get only clusters with minimum number of members
min_no_members <- 15 #TODO: change?

large_cluster_ids <- cluster_counts %>%
  filter(count >= min_no_members) %>% 
  .$cluster

# create an empty list to store the pie charts
pie_charts <- list()

# loop through each cluster group and create a pie chart
for (i in large_cluster_ids) {
  # filter the data for the current cluster group
  subset_df <- subset(max_groups_phys, cluster == i)
  
  # create a pie chart for the current cluster group
  pie_chart <- ggplot(subset_df, aes(x = "", fill = kidney_disease_group_short)) +
    geom_bar(width = 1) +
    coord_polar(theta = "y") +
    labs(title = paste("Cluster", i))
  
  # add the pie chart to the list
  pie_charts[[i]] <- pie_chart
}

# calculate the relative frequency of each cluster
cluster_freq <- prop.table(table(max_groups_phys$cluster))

# define a custom color palette for kidney_disease_group_short
custom_colors <- c("tubulopathy" = "Red", 
                   "glomerulopathy" = "Green", 
                   "cancer" = "Blue", 
                   "cakut" = "Purple", 
                   "cyst_cilio" = "Orange", 
                   "complement" = "Yellow",
                   "nephrocalcinosis" = "Grey") 

# create a data frame for pie chart creation
pie_data <- data.frame(cluster = names(cluster_freq), freq = cluster_freq)

# create and save a pie chart for each cluster group
for (i in large_cluster_ids) {
  subset_df <- subset(max_groups_phys, cluster == i)
  
  # calculate the total number of instances in the current cluster
  total_instances <- nrow(subset_df)
  
  pie_chart <- ggplot(subset_df, aes(x = "", fill = kidney_disease_group_short)) +
    geom_bar(width = 1) +
    scale_fill_manual(values = custom_colors) +
    coord_polar(theta = "y") +
    labs(title = paste("Physical network - Cluster", i, "- Total Instances:", total_instances))
  
  # define the filename based on cluster number
  filename <- paste0(output_dir, "/protein_interactions_phys_network_cluster_", i, "_pie_chart.", current_date, ".png")
  
  # save the pie chart as a PNG file
  ggsave(filename, plot = pie_chart, width = 5, height = 5)  
}
############################################

# TODO: get enrichment pathways etc
# group1 <- max_groups_full %>% filter(cluster == 1)
# 
# string_db_full$get_enrichment(group1$STRING_id) %>%
#   tibble() %>% View
# junk <- string_db_full$get_enrichment(group1$STRING_id) %>%
#   tibble()
