require(tidyverse)

#' Get TPM Expression Values from Descartes Atlas
#'
#' This function takes a character vector of gene names and a string representing
#' the tissue type as input, constructs the URL to fetch the CSV file from
#' Descartes Atlas, and returns the TPM expression values for the specified genes
#' in the given tissue as a tibble.
#'
#' @param genes A character vector representing the gene names.
#' @param tissue A character string representing the tissue type.
#'
#' @return A tibble with columns "Gene" and "TPM" containing the TPM expression 
#'         data for the specified genes in the given tissue.
#'
#' @examples
#' \dontrun{
#'   expression_data <- get_descartes_tpm_expression(c("IGF2", "MALAT1", "H19"), "kidney")
#'   print(expression_data)
#' }
#' \dontrun{
#'   expression_data <- get_descartes_tpm_expression(c("IGF2", "MALAT1", "H19"), "muscle")
#'   print(expression_data)
#' }
#'
#' @export
get_descartes_tpm_expression <- function(genes, tissue) {

  # Construct the URL based on the tissue input
  url <- paste0(
    "https://atlas.fredhutch.org/data/bbi/descartes/human_gtex/tables/tissue_tpm/",
    tissue,
    ".csv"
  )

  # Read the CSV file from the URL
  expression_data <- read_csv(url, col_names = c("Gene", "TPM"))

  # Filter the expression data to include only the input genes
  filtered_expression_data <- expression_data %>% 
    filter(Gene %in% genes)

  # Return the filtered expression data
  return(filtered_expression_data)
}


#' Get Percentage Of Cells That Express Genes from Descartes Atlas
#'
#' This function takes a character vector of gene names and a string representing
#' the tissue type as input, constructs the URL to fetch the CSV file from
#' Descartes Atlas, and returns the percentage of cells in the specified tissue
#' that express the given genes as a tibble.
#'
#' @param genes A character vector representing the gene names.
#' @param tissue A character string representing the tissue type.
#'
#' @return A tibble with columns "Gene" and "Percentage" containing the percentage
#'         of cells in the specified tissue that express the given genes.
#'
#' @examples
#' \dontrun{
#'   percentage_data <- get_descartes_cell_percentage_expression(
#'     c("MALAT1", "IGF2", "H19"),
#'     "kidney"
#'   )
#'   print(percentage_data)
#' }
#' \dontrun{
#'   percentage_data <- get_descartes_cell_percentage_expression(
#'     c("MALAT1", "IGF2", "H19"),
#'     "muscle"
#'   )
#'   print(percentage_data)
#' }
#'
#' @export
get_descartes_cell_percentage_expression <- function(genes, tissue) {

  # Construct the URL based on the tissue input
  url <- paste0(
    "https://atlas.fredhutch.org/data/bbi/descartes/human_gtex/tables/tissue_percentage/", 
    tissue, 
    ".csv"
  )

  # Read the CSV file from the URL
  percentage_data <- read_csv(url, col_names = c("Gene", "Percentage"))

  # Filter the percentage data to include only the input genes
  filtered_percentage_data <- percentage_data %>% 
    filter(Gene %in% genes)

  # Return the filtered percentage data
  return(filtered_percentage_data)
}
