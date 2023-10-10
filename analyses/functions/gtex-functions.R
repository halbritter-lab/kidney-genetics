require(httr)
require(jsonlite)
require(tidyverse)
require(purrr)
require(dplyr)
require(janitor)


#' Get Median Tissue Expression Levels from GTEx API
#'
#' This function takes a character vector of GENCODE gene identifiers as input,
#' makes a GET request to the GTEx API, and returns the median tissue
#' expression levels as a tibble.
#'
#' @param gencode_ids A character vector representing the GENCODE gene identifiers.
#' @param tissue_site_detail_ids A character vector of tissue site detail IDs to filter the output. Default is NULL.
#' @param wide_format A logical indicating whether to return the data in a wide format. Default is FALSE.
#'
#' @return A tibble with columns containing the median tissue expression data
#'         returned by the GTEx API. The format of the tibble (wide or long)
#'         can be controlled using the wide_format parameter.
#'
#' @examples
#' \dontrun{
#'   median_expression <- get_median_tissue_expression(c("ENSG00000008710.19", "ENSG00000118762.7"))
#'   print(median_expression)
#' }
#' \dontrun{
#'   median_expression <- get_median_tissue_expression(c("ENSG00000008710.19", "ENSG00000118762.7"), wide_format = TRUE)
#'   print(median_expression)
#' }
#'
#' @export
get_median_tissue_expression <- function(gencode_ids, tissue_site_detail_ids = NULL, wide_format = FALSE) {

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

  # Create a template dataframe with all gencodeIds from the input
  template_df <- tibble(gencodeId = gencode_ids)

  # Filter the data based on tissue_site_detail_ids
  if (!is.null(tissue_site_detail_ids)) {
    expression_tibble <- expression_tibble %>%
      filter(tissueSiteDetailId %in% tissue_site_detail_ids)
  }

  # Pivot to wide format if wide_format is TRUE
  if (wide_format) {
    expression_tibble <- expression_tibble %>%
      dplyr::select(gencodeId, median, tissueSiteDetailId) %>%
      pivot_wider(names_from = tissueSiteDetailId, values_from = median)
  }

  # Left join with the template to include all gencodeIds in the output
  output_tibble <- template_df %>%
    left_join(expression_tibble, by = "gencodeId")

  # Convert column names to snake case
  output_tibble <- janitor::clean_names(output_tibble)

  return(output_tibble)
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
#'        that can be queried in a single API request. Default is 5.
#' @param tissue_site_detail_ids A character vector of tissue site detail IDs to filter the output. Default is NULL.
#' @param wide_format A logical indicating whether to return the data in a wide format. Default is FALSE.
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
#'   median_expression <- get_multiple_median_tissue_expression(c("ENSG00000008710.19", "ENSG00000118762.7"), max_ids_per_request = 5)
#'   print(median_expression)
#' }
#'
#' @export
get_multiple_median_tissue_expression <- function(gencode_ids, max_ids_per_request = 5, tissue_site_detail_ids = NULL, wide_format = FALSE) {

  # Split the gencode_ids into chunks that fit within the API's limitations
  id_chunks <- split(gencode_ids, ceiling(seq_along(gencode_ids) / max_ids_per_request))

  # Use purrr::map to iterate over id_chunks and get a list of tibbles
  list_of_tibbles <- purrr::map(id_chunks, ~get_median_tissue_expression(.x, tissue_site_detail_ids = tissue_site_detail_ids, wide_format = wide_format))

  # Use purrr::reduce to bind all tibbles into a single tibble
  aggregated_data <- purrr::reduce(list_of_tibbles, dplyr::bind_rows)

  # Return the aggregated tibble
  return(aggregated_data)
}


#' Get GENCODE IDs from GTEx API
#'
#' This function takes a list of gene symbols as input, makes a GET request
#' to the GTEx API, and returns the first GENCODE ID encountered for each
#' gene symbol as a tibble. If the API returns an empty response, NA is returned
#' for the gencode_id.
#'
#' @param gene_symbols A character vector representing the gene symbols.
#'
#' @return A tibble with columns: gene_symbol and gencode_id.
#'
#' @examples
#' \dontrun{
#'   gencode_ids <- get_gencode_ids(c("PKD1", "PKD2"))
#'   print(gencode_ids)
#' }
#' @export
get_gencode_ids <- function(gene_symbols) {

  # Prepare the API URL
  base_url <- "https://gtexportal.org/api/v2/reference/gene?"
  ids_query <- paste0("geneId=", gene_symbols, collapse = "&")
  full_url <- paste0(base_url, ids_query)

  # Make the API request
  response <- GET(full_url)

  # Check if request was successful
  if (http_status(response)$category != "Success") {
    stop("Failed to retrieve data from GTEx API.")
  }

  # Parse the JSON response
  parsed_data <- fromJSON(content(response, "text", encoding = "UTF-8"))

  # Create a tibble from the input gene symbols
  output_tibble <- tibble(gene_symbol = gene_symbols)

  # Check if the response contains data
  if (length(parsed_data$data) > 0) {
    # Create a tibble from the parsed data
    response_tibble <- as_tibble(parsed_data$data)

    # Join the input tibble with the response tibble to ensure all input gene
    # symbols appear in the output and select only the first GENCODE ID for
    # each gene symbol if multiple are returned.
    output_tibble <- output_tibble %>%
      left_join(response_tibble, by = c("gene_symbol" = "geneSymbol")) %>%
      group_by(gene_symbol) %>%
      slice(1) %>%
      ungroup() %>%
      dplyr::select(gene_symbol, gencode_id = gencodeId)
  } else {
    # If the response is empty, add a column with NA values for the gencode_id
    output_tibble <- output_tibble %>%
      mutate(gencode_id = NA_character_)
  }

  return(output_tibble)
}


#' Fetch GENCODE IDs from GTEx API for Multiple Gene Symbols
#'
#' This function takes a character vector of gene symbols as input, splits it
#' into appropriate chunks to fit the API input size limitations, makes a GET
#' request for each chunk of gene symbols to the GTEx API, and returns the
#' GENCODE IDs as a tibble.
#'
#' @param gene_symbols A character vector representing the gene symbols.
#' @param max_symbols_per_request An integer specifying the maximum number of
#'        gene symbols that can be queried in a single API request. Default is 5.
#'
#' @return A tibble with rows for each gene and a column containing the
#'         GENCODE ID returned by the GTEx API.
#'
#' @examples
#' \dontrun{
#'   gencode_ids <- get_multiple_gencode_ids(c("PKD1", "PKD2"))
#'   print(gencode_ids)
#' }
#' \dontrun{
#'   gencode_ids <- get_multiple_gencode_ids(c("PKD1", "PKD2"), max_symbols_per_request = 5)
#'   print(gencode_ids)
#' }
#' @export
get_multiple_gencode_ids <- function(gene_symbols, max_symbols_per_request = 5) {

  # Split the gene_symbols into chunks that fit within the API's limitations
  symbol_chunks <- split(gene_symbols, ceiling(seq_along(gene_symbols) / max_symbols_per_request))

  # Use purrr::map to iterate over symbol_chunks and get a list of tibbles
  list_of_tibbles <- purrr::map(symbol_chunks, get_gencode_ids)

  # Use purrr::reduce to bind all tibbles into a single tibble
  aggregated_data <- purrr::reduce(list_of_tibbles, dplyr::bind_rows)

  # Return the aggregated tibble
  return(aggregated_data)
}
