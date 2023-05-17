#### This file holds analyses functions for PubTator requests


#' This function finds genes from a PubTator API request.
#'
#' @param query Character. The query string to be searched in PubTator.
#' @param page Numeric. The page number for the search results.
#' @param filter_type Character. The type of entity to be filtered from the
#' results. Default is "Gene".
#'
#' @return A tibble with the PubMed ID, text, text identifier, text part,
#' source, and text part count for each entry matching the filter type.
#' Returns NULL if no results found.
#'
#' @examples
#' genes <- pubtator_genes_in_request(query = "BRCA1", page = 1)
#'
#' @export
pubtator_genes_in_request <- function(query, page, filter_type = "Gene") {
  url <- paste0("https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/search?q=", query, "&page=", page)
  search_request <- fromJSON(URLencode(url), flatten = TRUE)

  # TODO: implement retry if error
  # TODO: also return NCBI gene ID

  search_results_tibble <- search_request$results %>%
    as_tibble()

  search_results_filtered <- search_results_tibble %>%
    {if (!("accessions" %in% colnames(.))) add_column(., accessions = "NULL") else .} %>%
    filter(accessions != "NULL")

  if (nrow(search_results_filtered) > 0)
  {
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