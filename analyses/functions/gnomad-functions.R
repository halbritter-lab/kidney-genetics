require(httr)
require(jsonlite)
require(tidyverse)
require(purrr)

#' Fetch Gene Data from gnomAD GraphQL API
#'
#' This function takes an Ensemble gene identifier as input, makes a POST request
#' to the gnomAD GraphQL API, and returns the gene data as a tibble.
#'
#' @param ensemble_id A character string representing the Ensemble gene identifier.
#'
#' @return A tibble with columns for each field returned by the API. This includes
#'         fields such as gene_id, gene_version, symbol, and others, as well as nested
#'         fields under gnomad_constraint.
#'
#' @examples
#' \dontrun{
#'   gene_data <- get_gene_data_from_gnomad("ENSG00000008710")
#'   print(gene_data)
#' }
#'
#' @export
get_gene_data_from_gnomad <- function(ensemble_id) {

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

  # Send POST request
  response <- POST(api_url, body = body, encode = "json")

  # Check if request was successful
  if (response$status_code == 200) {
    # Parse the JSON response
    data <- fromJSON(content(response, "text", encoding = "UTF-8"))

    # Convert the nested list to a flat data frame
    df <- as.data.frame(data$data$gene)

    # Convert the data frame to a tibble and return
    return(as_tibble(df))
  } else {
    stop("Request failed with status code ", response$status_code)
  }
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
#'   variants_data <- getClinVarVariants("ENSG00000008710")
#'   print(variants_data)
#' }
#'
#' @export
getClinVarVariants <- function(ensemble_id) {

  # API URL
  api_url <- "https://gnomad.broadinstitute.org/api/"

  # Request Body
  body <- list(
    query = paste0(
      '{
        gene(gene_id: "', ensemble_id,'", reference_genome: GRCh37) {
          clinvar_variants {
            clinical_significance
            clinvar_variation_id
            gold_stars
            hgvsc
            hgvsp
            in_gnomad
            major_consequence
            pos
            review_status
            transcript_id
            variant_id
          }
        }
      }'
    )
  )

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
  } else {
    stop("Request failed with status code ", response$status_code)
  }
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
