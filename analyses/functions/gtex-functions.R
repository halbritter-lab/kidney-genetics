require(httr)
require(jsonlite)
require(tidyverse)
require(purrr)


#' Get Median Tissue Expression Levels from GTEx API
#'
#' This function takes a list of GENCODE gene identifiers as input,
#' makes a GET request to the GTEx API, and returns the median tissue
#' expression levels as a tibble.
#'
#' @param gencode_ids A character vector representing the GENCODE gene identifiers.
#'
#' @return A tibble with columns containing the median tissue expression data
#'         returned by the GTEx API.
#'
#' @examples
#' \dontrun{
#'   median_expression <- get_median_tissue_expression("ENSG00000008710.19", "ENSG00000118762.7"))
#'   print(median_expression)
#' }
#'
#' @export
get_median_tissue_expression <- function(gencode_ids) {

  # Prepare the API URL
  base_url <- "https://gtexportal.org/api/v2/expression/medianGeneExpression?"
  ids_query <- paste0("gencodeId=", gencode_ids, collapse = "&")
  full_url <- paste0(base_url, ids_query)

  # Make the API request
  response <- GET(full_url)

  # Check if request was successful
  if (http_status(response)$category != "Success") {
    stop("Failed to retrieve data from GTEx API.")
  }

  # Parse the JSON response
  parsed_data <- fromJSON(content(response, "text", encoding = "UTF-8"))

  # Convert the parsed JSON data to a tibble
  expression_tibble <- as_tibble(parsed_data$data)

  return(expression_tibble)
}


#' Fetch Median Tissue Expression Levels from GTEx API for Multiple Genes
#'
#' This function takes a character vector of GENCODE gene identifiers as input,
#' splits it into appropriate chunks to fit the API input size limitations, 
#' makes a GET request for each chunk of gene identifiers to the GTEx API, 
#' and returns the median tissue expression levels as a tibble.
#'
#' @param gencode_ids A character vector representing the GENCODE gene identifiers.
#' @param max_ids_per_request An integer specifying the maximum number of identifiers
#'        that can be queried in a single API request. Default is 50.
#'
#' @return A tibble with rows for each gene and columns containing the median tissue
#'         expression data returned by the GTEx API.
#'
#' @examples
#' \dontrun{
#'   median_expression <- get_multiple_median_tissue_expression(c("ENSG00000008710.19", "ENSG00000118762.7"))
#'   print(median_expression)
#' }
#' \dontrun{
#'   median_expression <- get_multiple_median_tissue_expression(c()"ENSG00000008710.19", "ENSG00000118762.7"), max_ids_per_request = 25)
#'   print(median_expression)
#' }
#'
#' @export
get_multiple_median_tissue_expression <- function(gencode_ids, max_ids_per_request = 50) {

  # Split the gencode_ids into chunks that fit within the API's limitations
  id_chunks <- split(gencode_ids, ceiling(seq_along(gencode_ids) / max_ids_per_request))

  # Use purrr::map to iterate over id_chunks and get a list of tibbles
  list_of_tibbles <- purrr::map(id_chunks, get_median_tissue_expression)

  # Use purrr::reduce to bind all tibbles into a single tibble
  aggregated_data <- purrr::reduce(list_of_tibbles, dplyr::bind_rows)

  # Return the aggregated tibble
  return(aggregated_data)
}
