#### This file holds analyses functions for the ensembl database name using biomart

#' Retrieve gene coordinates in BED format from gene symbols
#'
#' This function retrieves the gene coordinates in BED format for the given gene
#' symbols. The coordinates are obtained from the specified reference genome.
#'
#' @param gene_symbols A vector or tibble containing the gene symbols.
#' @param reference The reference genome to use (default: "hg19").
#'
#' @return A tibble with the gene symbols and their corresponding coordinates in BED format.
#'
#' @examples
#' gene_symbols <- c("ARID1B ", "GRIn2B", "NAA10")
#' gene_coordinates_from_symbol(gene_symbols, reference = "hg19")
#'
#' @export
gene_coordinates_from_symbol <- function(gene_symbols, reference = "hg19") {
  gene_symbol_list <- as_tibble(gene_symbols) %>%
    dplyr::select(hgnc_symbol = value)

  # define mart
  mart_hg19 <- useMart("ensembl", host="grch37.ensembl.org")
  mart_hg19 <- useDataset("hsapiens_gene_ensembl", mart_hg19)

  mart_hg38 <- useMart("ensembl", host="ensembl.org")
  mart_hg38 <- useDataset("hsapiens_gene_ensembl", mart_hg38)

  if (reference == "hg19") {
    mart <- useMart("ensembl", host = "grch37.ensembl.org")
    mart <- useDataset("hsapiens_gene_ensembl", mart_hg19)
  } else {
    mart <- useMart("ensembl", host = "ensembl.org")
    mart <- useDataset("hsapiens_gene_ensembl", mart_hg38)
  }

  attributes <- c("hgnc_symbol", "chromosome_name", "start_position", "end_position")
  filters <- c("hgnc_symbol")

  values <- list(hgnc_symbol = gene_symbol_list$hgnc_symbol)

  gene_coordinates_hg19 <- getBM(attributes=attributes, filters=filters, values=values, mart=mart) %>%
    group_by(hgnc_symbol) %>%
    summarise(hgnc_symbol = max(hgnc_symbol), chromosome_name = max(chromosome_name), start_position = max(start_position), end_position = max(end_position)) %>%
    mutate(bed_format = paste0("chr", chromosome_name, ":", start_position, "-", end_position)) %>%
    dplyr::select(hgnc_symbol, bed_format)

  gene_symbol_list_return <- gene_symbol_list %>%
  left_join(gene_coordinates_hg19, by = ("hgnc_symbol"))

  return(gene_symbol_list_return)
}


#' Retrieve gene coordinates in BED format from Ensembl IDs
#'
#' This function retrieves the gene coordinates in BED format for the given Ensembl
#' gene IDs. The coordinates are obtained from the specified reference genome.
#'
#' @param ensembl_id A vector or tibble containing the Ensembl gene IDs.
#' @param reference The reference genome to use (default: "hg19").
#'
#' @return A tibble with the Ensembl gene IDs and their corresponding coordinates in BED format.
#'
#' @examples
#' ensembl_id <- c("ENSG00000123456", "ENSG00000123457", "ENSG00000123458")
#' gene_coordinates_from_ensembl(ensembl_id, reference = "hg19")
#'
#' @export
gene_coordinates_from_ensembl <- function(ensembl_id, reference = "hg19") {
  ensembl_id_list <- as_tibble(ensembl_id) %>%
    dplyr::select(ensembl_gene_id = value)

  # define mart
  mart_hg19 <- useMart("ensembl", host="grch37.ensembl.org")
  mart_hg19 <- useDataset("hsapiens_gene_ensembl", mart_hg19)

  mart_hg38 <- useMart("ensembl", host="ensembl.org")
  mart_hg38 <- useDataset("hsapiens_gene_ensembl", mart_hg38)

  if (reference == "hg19") {
    mart <- useMart("ensembl", host = "grch37.ensembl.org")
    mart <- useDataset("hsapiens_gene_ensembl", mart_hg19)
  } else {
    mart <- useMart("ensembl", host = "ensembl.org")
    mart <- useDataset("hsapiens_gene_ensembl", mart_hg38)
  }

  attributes <- c("ensembl_gene_id", "chromosome_name", "start_position", "end_position")
  filters <- c("ensembl_gene_id")

  values <- list(ensembl_gene_id = ensembl_id_list$ensembl_gene_id)

  gene_coordinates_hg19 <- getBM(attributes=attributes, filters=filters, values=values, mart=mart) %>%
    group_by(ensembl_gene_id) %>%
    summarise(ensembl_gene_id = max(ensembl_gene_id), chromosome_name = max(chromosome_name), start_position = max(start_position), end_position = max(end_position)) %>%
    mutate(bed_format = paste0("chr", chromosome_name, ":", start_position, "-", end_position)) %>%
    dplyr::select(ensembl_gene_id, bed_format)

  ensembl_id_list_return <- ensembl_id_list %>%
  left_join(gene_coordinates_hg19, by = ("ensembl_gene_id"))

  return(ensembl_id_list_return)
}