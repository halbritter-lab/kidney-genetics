require(tidyverse)


#' Compute Protein Interactions Using STRINGdb Data
#'
#' This function takes a character vector of protein identifiers in the format
#' of STRINGdb, downloads and caches the STRINGdb protein interaction data,
#' computes interaction scores, and returns the results as a tibble.
#'
#' @param protein_ids A character vector representing the protein identifiers
#'        in STRINGdb format.
#' @param refresh A logical indicating whether to refresh the cached STRINGdb
#'        file. Default is FALSE.
#' @param max_age_months An integer specifying the maximum age of the cached
#'        file in months. If the file is older, it will be re-downloaded.
#'        Default is 1.
#' @param directory A string representing the directory where the STRINGdb file
#'        is cached. Default is "./".
#'
#' @return A tibble with rows for each protein and columns containing the summed
#'         and normalized interaction scores, and a string representation of
#'         interacting proteins and scores.
#'
#' @examples
#' \dontrun{
#'   protein_ids <- c("9606.ENSP00000000233", "9606.ENSP00000257770", 
#'                    "9606.ENSP00000226004")
#'   interaction_results <- compute_protein_interactions(protein_ids)
#'   print(interaction_results)
#' }
#'
#' @export
compute_protein_interactions <- function(protein_ids, refresh = FALSE, 
                                         max_age_months = 1, directory = "./") {

  # Define the URL and filename prefix for the STRINGdb file
  stringdb_url <- paste0(
    "https://stringdb-downloads.org/download/protein.physical.links.v12.0/",
    "9606.protein.physical.links.v12.0.txt.gz"
  )
  filename_prefix <- "stringdb_protein_links"

  # Ensure the directory ends with a slash
  if (substr(directory, nchar(directory), nchar(directory)) != "/") {
    directory <- paste0(directory, "/")
  }

  # Check if the file needs to be downloaded
  if (refresh || !check_file_age(filename_prefix, directory, max_age_months)) {
    # Define the filename with the current date
    filename <- paste0(directory, filename_prefix, ".", Sys.Date(), ".txt.gz")

    # Download the file
    download_status <- try(download.file(stringdb_url, filename, mode = "wb"), 
                           silent = TRUE)
    if (inherits(download_status, "try-error")) {
      stop("Failed to download the STRINGdb file.")
    }
  } else {
    # Use the newest cached file
    filename <- get_newest_file(filename_prefix, directory)
  }

  # Read the data from the file
  interaction_data <- read_table(filename, col_types = "ccn")

  # Filter the data to include only interactions between the input proteins
  filtered_data <- interaction_data %>% 
    filter((protein1 %in% protein_ids & protein2 %in% protein_ids) |
             (protein2 %in% protein_ids & protein1 %in% protein_ids))

  # Sum the scores for each protein
  summed_scores <- filtered_data %>% 
    group_by(protein1) %>% 
    summarise(sum_score = sum(combined_score)) %>% 
    ungroup()

  # Percentile normalize the summed scores
  summed_scores$normalized_score <- percent_rank(summed_scores$sum_score)

  # Create a string of interacting proteins and scores for each protein1
  interaction_strings <- filtered_data %>% 
    group_by(protein1) %>%
    summarise(interaction_string = paste(paste(protein2, combined_score, 
                                                sep = ":"), 
                                          collapse = ", ")) %>%
    ungroup()

  # Create a base tibble with all input protein IDs
  base_tibble <- tibble(protein1 = protein_ids)

  # Join interaction_strings and summed_scores with base_tibble
  result_tibble <- base_tibble %>%
    left_join(summed_scores, by = "protein1") %>%
    left_join(interaction_strings, by = "protein1")

  # Return the results
  return(result_tibble)
}
