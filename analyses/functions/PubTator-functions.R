#### This file holds analyses functions for PubTator requests


#' Retrieve Gene-Related Data from PubTator API v2
#'
#' This function queries the PubTator v2 API to retrieve gene-related data based on a given query string.
#' It constructs the API request URL using parameters for the base URL, endpoint, and query formatting.
#' The function handles pagination and can retry requests in case of temporary API issues.
#'
#' @param query Character: The search query string for PubTator.
#' @param page Numeric: The page number for the API response (for pagination).
#' @param max_retries Numeric: Maximum number of retries for the API request in case of failure. Defaults to 3.
#' @param api_base_url Character: Base URL of the PubTator API. Defaults to the v2 API URL.
#' @param endpoint Character: API endpoint for the search query. Defaults to "search".
#' @param query_parameter Character: URL parameter for the search query. Defaults to "?q=".
#' @param filter_type Character: Type of entity to filter from the results. Defaults to "Gene".
#'
#' @return A tibble containing gene-related data from PubTator. Each row represents a unique gene-related entry
#'   with fields including PubMed ID (pmid), text of the gene mention (text), gene identifier (text_identifier),
#'   part of the article where the gene is mentioned (text_part), source of the mention (source), and count of 
#'   text parts mentioning the gene (text_part_count). Returns NULL if no results are found.
#'
#' @examples
#' genes <- pubtator_v2_genes_in_request(query = "BRCA1", page = 10)
#' # Example output:
#' # A tibble: 93 × 6
#' #      pmid text     text_identifier text_part source text_part_count
#' #     <int> <chr>    <chr>           <chr>     <chr>            <int>
#' # 1 37794204 BRCA     672             title     Gene                 1
#' # 2 37794204 BRCA1/2  672;675         title     Gene                 1
#' # ... with more rows
#'
#' @export
pubtator_genes_in_request <- function(query, page, filter_type = "Gene", max_retries = 3) {
  # TODO: also return NCBI gene ID

  # define URL
  url <- paste0("https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/search?q=", query, "&page=", page)
  retries <- 0

  while(retries <= max_retries) {
    # Attempt to get data
    tryCatch({
      search_request <- fromJSON(URLencode(url), flatten = TRUE)

      # If successful, exit loop
      break

    }, error = function(e) {
      # If there's an error, increment the retry count and print a warning
      retries <- retries + 1
      if(retries <= max_retries) {
        warning(paste("Attempt", retries, "failed. Retrying..."))
      }
    })
  }

  # If max retries exhausted and still failed, throw an error
  if(retries > max_retries) {
    stop("Failed to fetch data after", max_retries, "attempts.")
  }

  # Continue processing the data as before
  search_results_tibble <- search_request$results %>%
    as_tibble()

  search_results_filtered <- search_results_tibble %>%
    {if (!("accessions" %in% colnames(.))) add_column(., accessions = "NULL") else .} %>%
    filter(accessions != "NULL")

  if (nrow(search_results_filtered) > 0) {
    search_results <- search_results_filtered %>%
      select(pmid, passages) %>%
      unnest(passages) %>%
      select(pmid, text_part = infons.type, annotations) %>%
      rowwise() %>%
      mutate(empty = is_empty(annotations)) %>%
      ungroup() %>%
      filter(annotations != "NULL" & !empty) %>%
      unnest(annotations, names_repair = "universal") %>%
      select(pmid, text_part, text, type = infons.type, text_identifier = infons.identifier) %>%
      unique() %>%
      filter(type == filter_type) %>%
      group_by(pmid, text, text_identifier) %>%
      summarise(text_part = paste(unique(text_part), collapse = " | "),
        source = paste(unique(type), collapse = " | "),
        text_part_count = n(),
        .groups = "keep") %>%
      ungroup()

    return(search_results)
  } else {
    return(NULL)
  }
}


#' This function returns the number of pages for a request to the NCBI PubTator
#' API.
#'
#' @param query Character. The query string to be searched in PubTator.
#'
#' @return Numeric. The total number of pages available for the given query.
#'
#' @examples
#' pages <- pubtator_pages_request(query = "BRCA1")
#'
#' @export
pubtator_pages_request <- function(query) {
  url <- paste0("https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/search?q=", query, "&page=", 1)
  search_request <- fromJSON(URLencode(url), flatten = TRUE)

  return(search_request$total_pages)
}