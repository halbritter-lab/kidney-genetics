#### This file holds analyses functions for NCBI datasets v2 API request


#' This function uses the NCBI datasets v2 API to return gene information based
#' on a NCBI gene ID.
#'
#' @param input A vector, list or other object that can be coerced to a tibble
#' containing NCBI gene IDs.
#' @param api_key Character. Your NCBI API key. By default it uses the value of
#' NCBI_API_KEY.
#' @param request_max Numeric. The maximum number of requests to be made at a
#' time to the API. Default is 20.
#'
#' @return A tibble with the input NCBI gene IDs and associated gene information.
#'
#' @examples
#' gene_info <- gene_info_from_gene_id(input = c("GeneID1", "GeneID2"),
#' api_key = "your_api_key")
#'
#' @export
gene_info_from_gene_id <- function(input, api_key = NCBI_API_KEY, request_max = 20) {
  ep_base_url <- "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/gene/id/"
  ep_fields_args <- "?table_fields=gene-id&table_fields=gene-type&table_fields=description"

  # TODO: implement retry if empty
  # TODO: implement error handling if input not a number or NA

  input_list <- as_tibble(input, .name_repair = "unique") %>%
    mutate(value = as.character(value))
  input_list_unique <- as_tibble(input, .name_repair = "unique") %>%
    unique()

  row_number <- nrow(input_list_unique)
  groups_number <- ceiling(row_number/request_max)

  input_list_result <- input_list_unique %>%
    mutate(group = sample(1:groups_number, row_number, replace=T)) %>%
    group_by(group) %>%
    summarise(query = paste(unique(value), collapse = ",")) %>%
    ungroup() %>%
    mutate(url = paste0(ep_base_url,
      query,
      ep_fields_args)) %>%
    rowwise() %>%
    mutate(results = list(fromJSON(content(GET(url, add_headers("api-key" = api_key), accept_json()), "text", encoding = "UTF-8"))$reports$gene %>% 
      as_tibble() %>%
      unnest(annotations) %>%
      unnest(nomenclature_authority) %>%
      unnest(ensembl_gene_ids, keep_empty = TRUE) %>%
      select(gene_id, symbol, tax_id, taxname, identifier, ensembl_gene_ids) %>%
      unique())) %>%
    select(-group, - query, -url) %>%
    unnest(results)

  input_list_return <- input_list %>%
    left_join(input_list_result, by = c("value" = "gene_id"))

  return(input_list_return)
}


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
    {if( !("accessions" %in% colnames(.))) add_column(., accessions = "NULL") else .} %>%
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