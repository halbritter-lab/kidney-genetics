require(httr)
require(jsonlite)
require(tidyverse)
require(purrr)

get_gene_data_from_gnomad <- function(ensemble_id) {

  # Recursive function to replace NULL with NA
  replace_null_with_na <- function(x) {
    if (is.list(x)) {
      lapply(x, replace_null_with_na)
    } else if (is.null(x)) {
      NA
    } else {
      x
    }
  }

  # API URL
  api_url <- "https://gnomad.broadinstitute.org/api/"

  # Request Body
  body <- list(
    query = paste0(
      '{
        gene(gene_id: "', ensemble_id, '", reference_genome: GRCh37) {
          reference_genome
          gene_id
          gene_version
          symbol
          gencode_symbol
          name
          canonical_transcript_id
          hgnc_id
          ncbi_id
          omim_id
          chrom
          start
          stop
          strand
          gnomad_constraint {
            exp_lof
            exp_mis
            exp_syn
            obs_lof
            obs_mis
            obs_syn
            oe_lof
            oe_lof_lower
            oe_lof_upper
            oe_mis
            oe_mis_lower
            oe_mis_upper
            oe_syn
            oe_syn_lower
            oe_syn_upper
            lof_z
            mis_z
            syn_z
            pLI
          }
        }
      }'
    )
  )

  # Set maximum number of retries
  max_retries <- 10

  # Initial waiting time in seconds
  wait_time <- 0.1

  for (i in 1:max_retries) {
    # Send POST request
    response <- POST(api_url, body = body, encode = "json")

    # Check if request was successful
    if (response$status_code == 200) {
      # Parse the JSON response
      data <- fromJSON(content(response, "text", encoding = "UTF-8"))

      # Replace NULL with NA
      data <- replace_null_with_na(data)

      # Convert the nested list to a flat data frame
      df <- as.data.frame(data$data$gene)

      # Convert the data frame to a tibble and return
      return(as_tibble(df))
    } else if (response$status_code == 429) {
      # Wait and then continue to next iteration (retry)
      Sys.sleep(wait_time)

      # Exponential backoff
      wait_time <- wait_time * 2 * i
    } else {
      stop("Request failed with status code ", response$status_code)
    }
  }

  # If code execution reached this point, all retries have failed
  stop("Request failed after ", max_retries, " retries")
}


#' Fetch ClinVar Variants Data from gnomAD GraphQL API
#'
#' This function takes an Ensemble gene identifier as input, makes a POST request
#' to the gnomAD GraphQL API, and returns the ClinVar variants data as a tibble.
#'
#' @param ensemble_id A character string representing the Ensemble gene identifier.
#'
#' @return A tibble with columns for each field returned by the API. This includes
#'         fields such as clinical_significance, clinvar_variation_id, and others, as well as nested
#'         fields under gnomad and exome/genome.
#'
#' @examples
#' \dontrun{
#'   variants_data <- get_clinvar_variants("ENSG00000008710")
#'   print(variants_data)
#' }
#'
#' @export
get_clinvar_variants <- function(ensemble_id) {

  # API URL
  api_url <- "https://gnomad.broadinstitute.org/api/"

  # Request Body
  body <- list(
    query = paste0(
      '{
        gene(gene_id: "', ensemble_id,'", reference_genome: GRCh37) {
          clinvar_variants {
            transcript_id
            hgvsc
            hgvsp
            major_consequence
            pos
            clinical_significance
            clinvar_variation_id
            variant_id
            gold_stars
            review_status
            in_gnomad
          }
        }
      }'
    )
  )

  # Set maximum number of retries
  max_retries <- 10

  # Initial waiting time in seconds
  wait_time <- 0.1

  for (i in 1:max_retries) {
    # Send POST request
    response <- POST(api_url, body = body, encode = "json")

    # Check if request was successful
    if (response$status_code == 200) {
      # Parse the JSON response
      data <- fromJSON(content(response, "text", encoding = "UTF-8"), flatten = TRUE)

      # Convert the nested list to a flat data frame
      df <- as.data.frame(data$data$gene)

      # Convert the data frame to a tibble and return
      return(as_tibble(df))
    } else if (response$status_code == 429) {
      # Wait and then continue to next iteration (retry)
      Sys.sleep(wait_time)

      # Exponential backoff
      wait_time <- wait_time * 2 * i
    } else {
      stop("Request failed with status code ", response$status_code)
    }
  }

  # If code execution reached this point, all retries have failed
  stop("Request failed after ", max_retries, " retries")
}


#' Fetch Gene Data from gnomAD GraphQL API for multiple genes
#'
#' This function takes a vector of Ensemble gene identifiers as input,
#' makes a POST request for each gene identifier
#' to the gnomAD GraphQL API, and returns the gene data as a tibble.
#'
#' @param ensemble_ids A character vector representing the Ensemble gene identifiers.
#'
#' @return A tibble with rows for each gene and columns for each field returned by the API.
#'         This includes fields such as gene_id, gene_version, symbol, and others, 
#'         as well as nested fields under gnomad_constraint.
#'
#' @examples
#' \dontrun{
#'   gene_data <- get_multiple_gene_data_from_gnomad(c("ENSG00000008710", "ENSG00000012048"))
#'   print(gene_data)
#' }
#'
#' @export
get_multiple_gene_data_from_gnomad <- function(ensemble_ids) {
  # Use purrr::map to iterate over ensemble_ids and get a list of tibbles
  list_of_tibbles <- purrr::map(ensemble_ids, get_gene_data_from_gnomad)

  # Use purrr::reduce to bind all tibbles into a single tibble
  gene_data <- purrr::reduce(list_of_tibbles, dplyr::bind_rows)

  # Return the combined tibble
  return(gene_data)
}


#' Fetch ClinVar Variants Data from gnomAD GraphQL API for multiple genes
#'
#' This function takes a vector of Ensemble gene identifiers as input,
#' makes a POST request for each gene identifier
#' to the gnomAD GraphQL API, and returns a list of ClinVar variants data as tibbles.
#'
#' @param ensemble_ids A character vector representing the Ensemble gene identifiers.
#'
#' @return A list with each element being a tibble with rows for each variant and columns for each field returned by the API. 
#'         This includes fields such as clinical_significance, clinvar_variation_id, 
#'         and others, as well as nested fields under gnomad and exome/genome.
#'
#' @examples
#' \dontrun{
#'   variants_data <- get_multiple_clinvar_variants(c("ENSG00000008710", "ENSG00000012048"))
#'   print(variants_data)
#' }
#'
#' @export
get_multiple_clinvar_variants <- function(ensemble_ids) {
  # Use purrr::map to iterate over ensemble_ids and get a list of tibbles
  list_of_tibbles <- purrr::map(ensemble_ids, get_clinvar_variants)

  # Return the list of tibbles
  return(list_of_tibbles)
}

