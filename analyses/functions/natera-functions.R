#### This file holds analyses functions for querying the natera web page

#' Get Nonce from Renasight Genetic Testing Page
#'
#' This function opens a curl handle, fetches the content of the Renasight Genetic Testing 
#' webpage, and extracts a nonce (a "number used once") using a specific xpath.
#' 
#' @param input_url A character string. The URL to fetch the nonce from.
#' Default: "https://www.natera.com/organ-health/renasight-genetic-testing/gene-conditions-list/"
#'
#' @return A character string. The extracted nonce from the Renasight Genetic Testing webpage.
#' 
#' @examples
#' natera_renasight_get_nonce()
#'
#' @export
natera_renasight_get_nonce <- function(input_url = "https://www.natera.com/organ-health/renasight-genetic-testing/gene-conditions-list/") {
  # open curl handle
  h <- new_handle()

  # fetch page
  r <- curl_fetch_memory(input_url)

  # find nonce using xpath
  nonce <- rawToChar(r$content) %>%
    read_html() %>%
    html_nodes(xpath = '//*[@id="_genescreening"]/@value') %>%
    html_text()

  # nonce
  return(nonce)
}


#' Get Last Page Number from Renasight Genetic Testing Page
#'
#' This function opens a curl handle, sets a POST request with necessary form data,
#' fetches the content of a webpage, and extracts the last page number using a specific xpath.
#'
#' @param input_url A character string. The URL to fetch the last page number from.
#' Default: "https://www.natera.com/wp-admin/admin-ajax.php"
#'
#' @return A character string. The extracted last page number from the webpage.
#' 
#' @examples
#' natera_renasight_get_last_page_number()
#'
#' @export
natera_renasight_get_last_page_number <- function(input_url = "https://www.natera.com/wp-admin/admin-ajax.php") {
  # get nonce
  natera_nonce <- natera_renasight_get_nonce()

  # open curl handle
  h <- new_handle()
    handle_setopt(h, customrequest = "POST")

  # set curl form
  handle_setform(h, .list = list(
    `action` = 'gene_screening_options',
    `page` = '1',
    `nonce` = natera_nonce
  ))

  # fetch page
  r <- curl_fetch_memory(input_url, h)

  # find last page using xpath
  last_page <- rawToChar(r$content) %>%
    read_html() %>%
    html_nodes(xpath = '//li[contains(@aria-label, "Last page")]/@data-page') %>%
    html_text()

  # return last page
  return(last_page)
}


#' Get Genes from Specific Page on Renasight Genetic Testing Website
#'
#' This function opens a curl handle, sets a POST request with necessary form data,
#' fetches the content of a specified webpage, and extracts a list of genes using a 
#' specific xpath. The function can also save the fetched page if a valid path is provided.
#' 
#' @param page_requested A numeric value or a character string. The specific page number
#' to fetch genes from.
#' @param input_url A character string. The URL to fetch genes from. 
#' Default: "https://www.natera.com/wp-admin/admin-ajax.php"
#' @param panel_name A character string. The name of the gene panel.
#' Default: "natera_renasight_comprehensive_kidney_gene_panel"
#' @param save_path A character string. The path to save fetched page. If the path doesn't
#' exist, the fetched page will not be saved. Default: "data/downloads"
#'
#' @return A character vector. The extracted genes from the specified webpage.
#' 
#' @examples
#' natera_renasight_get_genes_from_page(1)
#'
#' @export
natera_renasight_get_genes_from_page <- function(page_requested,
    input_url = "https://www.natera.com/wp-admin/admin-ajax.php",
    panel_name = "natera_renasight_comprehensive_kidney_gene_panel",
    save_path = "data/downloads") {
  # secure page_requested is a chara
  page_requested <- as.character(page_requested)

  # get nonce
  natera_nonce <- natera_renasight_get_nonce()

  # open curl handle
  h <- new_handle()
    handle_setopt(h, customrequest = "POST")

  # set curl form
  handle_setform(h, .list = list(
    `action` = 'gene_screening_options',
    `page` = page_requested,
    `nonce` = natera_nonce
  ))

  # fetch and save page if path exists
  if (file.exists(save_path)) {
    # generate file name
    creation_date <- strftime(as.POSIXlt(Sys.time(), "UTC", "%Y-%m-%dT%H:%M:%S"), "%Y-%m-%d")
    output_filename <- paste0(save_path, "/", panel_name, ".page_", page_requested, ".", creation_date, ".html")

    # make fetch disk request
    r <- curl_fetch_disk(input_url, output_filename, h)

    # assign content
    page_content <- r$content
  } else {  
    # make fetch memmory request
    r <- curl_fetch_memory(input_url, h)

    # assign content
    page_content <- rawToChar(r$content)
  }

  # find genes using xpath
  genes_in_page <- page_content %>%
    read_html() %>%
    html_nodes(xpath = '//div[contains(@class,"title-row__wrapper")]//span') %>%
    html_text()

  # return genes
  return(genes_in_page)
}