#### This file holds analyses functions for querying the blueprintgenetics sub panel pages

#' Extract Genes from BlueprintGenetics Panel Page
#'
#' This function extracts a list of genes from a specified BlueprintGenetics panel 
#' page. It finds the genes by identifying a specific HTML table element on the page.
#'
#' @param base_url A character string. The URL of the BlueprintGenetics panel page 
#' from which to extract genes.
#'
#' @return A character vector. The extracted genes from the specified BlueprintGenetics 
#' panel page.
#' 
#' @examples
#' blueprintgenetics_panel_extract_genes("https://www.blueprintgenetics.com/panels/panel-url/")
#'
#' @export
blueprintgenetics_panel_extract_genes <- function(base_url) {
	panel_page <- read_html(base_url)

	panel_genes <- (panel_page %>%
		html_nodes(xpath = '//table[contains(@class,"table mb-5")]') %>%
		html_table())[[1]] %>%
		select(Genes = Gene)

	return(panel_genes$Genes)
}