# load libraries
library(STRINGdb)
library(ggplot2)
# library(gridExtra)

# TODO: set in config?
string_db_version <- "12.0"
string_db_file_path <- "."

# create a directory to save the plots of the analysis
# TODO: change output_dir 
output_dir <- "plots"
dir.create(output_dir, showWarnings = TRUE)

# define a function to return the latest file with a given prefix in a given folder
get_latest_file <- function(folder_path, prefix) {
  # list files in the specified folder
  pattern1 <- paste0(prefix, "\\.\\d{4}-\\d{2}-\\d{2}\\.csv\\.gz")
  
  files <- list.files(folder_path, pattern = pattern1, full.names = TRUE)
  
  if (length(files) == 0) {
    return(NULL)
  }
  
  # extract dates from filenames
  pattern2 <- paste0(".*", prefix, "\\.(\\d{4}-\\d{2}-\\d{2})\\.csv\\.gz")
  dates <- gsub(pattern2, "\\1", files)
  dates <- as.Date(dates)
  
  # find the latest date
  latest_date <- max(dates)
  
  # convert the latest date back to the filename format
  latest_date_str <- format(latest_date, "%Y-%m-%d")
  latest_file <- files[dates == latest_date_str]
  
  return(latest_file)
}

# load annotated HGNC table
hgnc_annotated_path <- get_latest_file(folder_path = "/Users/nrank/Desktop/BioInf/halbritter/kidney-genetics/analyses/B_AnnotationHGNC/results",
                                       prefix = "non_alt_loci_set_coordinates")

hgnc_annotated <- read.csv(gzfile(hgnc_annotated_path)) %>% 
  dplyr::select(hgnc_id, symbol, STRING_id) %>% 
  mutate(hgnc_id = as.numeric(sub("HGNC:", "", hgnc_id)))

# load annotated high evidence genes
high_evidence_annotated_path <- get_latest_file(folder_path = "/Users/nrank/Desktop/BioInf/halbritter/kidney-genetics/analyses/C_AnnotateMergedTable/results",
                                                prefix = "C_high_evidence_annotated_csv_table") # TODO: change folder_path!!!

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

##### Analysis with FULL network ####
# download STRING db protein links 
string_url <- paste0("https://stringdb-downloads.org/download/protein.links.v", string_db_version, "/9606.protein.links.v", string_db_version, ".txt.gz")
download.file(string_url,
              destfile = paste0("9606.protein.links.v", string_db_version, ".txt.gz") ) 
# TODO: manual download required??

# instantiate new STRING db reference class 
string_db_full <- STRINGdb::STRINGdb$new(version = string_db_version,
                                    species = 9606,
                                    score_threshold = 100,
                                    input_directory = ".")

# get clusters
clusters_list_full <- string_db_full$get_clusters(unique(kid_groups$STRING_id),
                                        algorithm = "walktrap") # NOTE: some STRING_ids get lost here => no group for them???

clusters_tibble_full <- tibble(clusters_list_full) %>% 
  select(STRING_id = clusters_list_full) %>% 
  mutate(cluster = row_number()) %>% 
  unnest_longer(col = "STRING_id")

max_groups <- kid_groups %>%
  group_by(approved_symbol) %>%
  filter(hpo_id_group_p == max(hpo_id_group_p)) %>% 
  left_join(clusters_tibble, by = "STRING_id") %>% 
  filter(!is.na(cluster))
# NOTE: if genes have more than one group with same and highest hpo_id_group_p => keep all groups 

# group the data by 'cluster' and calculate the counts for each group
cluster_counts <- table(max_groups$cluster)

# create an empty list to store the pie charts
pie_charts <- list()

# loop through each cluster group and create a pie chart
for (i in unique(max_groups$cluster)) {
  # filter the data for the current cluster group
  subset_df <- subset(max_groups, cluster == i)
  
  # create a pie chart for the current cluster group
  pie_chart <- ggplot(subset_df, aes(x = "", fill = kidney_disease_group_short)) +
    geom_bar(width = 1) +
    coord_polar(theta = "y") +
    labs(title = paste("Cluster", i))
  
  # add the pie chart to the list
  pie_charts[[i]] <- pie_chart
}

# calculate the relative frequency of each cluster
cluster_freq <- prop.table(table(max_groups$cluster))

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
  subset_df <- subset(max_groups, cluster == pie_data$cluster[i])
  
  # calculate the total number of instances in the current cluster
  total_instances <- nrow(subset_df)
  
  pie_chart <- ggplot(subset_df, aes(x = "", fill = kidney_disease_group_short)) +
    geom_bar(width = 1) +
    scale_fill_manual(values = custom_colors) +
    coord_polar(theta = "y") +
    labs(title = paste("Full network - Cluster", pie_data$cluster[i], "- Total Instances:", total_instances))
  
  # define the filename based on cluster number
  filename <- paste0(output_dir, "/protein_interactions_full_network_cluster_", pie_data$cluster[i], "_pie_chart.png")
  
  # save the pie chart as a PNG file
  ggsave(filename, plot = pie_chart, width = 5, height = 5)  
}


##### Analysis with PHYSICAL network ####
# download STRING db protein links 
string_url <- paste0("https://stringdb-static.org/download/protein.physical.links.v", string_db_version, "/9606.protein.physical.links.v", string_db_version, ".txt.gz")
download.file(string_url,
              destfile = paste0("9606.protein.physical.links.v", string_db_version, ".txt.gz") ) 
# TODO: manual download required??

# instantiate new STRING db reference class 
string_db_phys <- STRINGdb::STRINGdb$new(version = string_db_version,
                                         species = 9606,
                                         score_threshold = 100,
                                         network_type = "physical",
                                         input_directory = ".")

# get clusters
clusters_list_phys <- string_db_phys$get_clusters(unique(kid_groups$STRING_id),
                                                  algorithm = "walktrap") # NOTE: some STRING_ids get lost here => no group for them???

clusters_tibble_phys <- tibble(clusters_list_phys) %>% 
  select(STRING_id = clusters_list_phys) %>% 
  mutate(cluster = row_number()) %>% 
  unnest_longer(col = "STRING_id")

max_groups <- kid_groups %>%
  group_by(approved_symbol) %>%
  filter(hpo_id_group_p == max(hpo_id_group_p)) %>% 
  left_join(clusters_tibble, by = "STRING_id") %>% 
  filter(!is.na(cluster))
# NOTE: if genes have more than one group with same and highest hpo_id_group_p => keep all groups 

# group the data by 'cluster' and calculate the counts for each group
cluster_counts <- table(max_groups$cluster)

# create an empty list to store the pie charts
pie_charts <- list()

# loop through each cluster group and create a pie chart
for (i in unique(max_groups$cluster)) {
  # filter the data for the current cluster group
  subset_df <- subset(max_groups, cluster == i)
  
  # create a pie chart for the current cluster group
  pie_chart <- ggplot(subset_df, aes(x = "", fill = kidney_disease_group_short)) +
    geom_bar(width = 1) +
    coord_polar(theta = "y") +
    labs(title = paste("Cluster", i))
  
  # add the pie chart to the list
  pie_charts[[i]] <- pie_chart
}

# calculate the relative frequency of each cluster
cluster_freq <- prop.table(table(max_groups$cluster))

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
  subset_df <- subset(max_groups, cluster == pie_data$cluster[i])
  
  # calculate the total number of instances in the current cluster
  total_instances <- nrow(subset_df)
  
  pie_chart <- ggplot(subset_df, aes(x = "", fill = kidney_disease_group_short)) +
    geom_bar(width = 1) +
    scale_fill_manual(values = custom_colors) +
    coord_polar(theta = "y") +
    labs(title = paste("Full network - Cluster", pie_data$cluster[i], "- Total Instances:", total_instances))
  
  # define the filename based on cluster number
  filename <- paste0(output_dir, "/protein_interactions_phys_network_cluster_", pie_data$cluster[i], "_pie_chart.png")
  
  # save the pie chart as a PNG file
  ggsave(filename, plot = pie_chart, width = 5, height = 5)  
}
