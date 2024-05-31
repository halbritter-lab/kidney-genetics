require(httr)
require(jsonlite)

#' Find the Last Page of PanelApp API Endpoint
#'
#' This function iteratively queries the PanelApp API endpoint, increasing the
#' page number until an "Invalid page" response is received. Once this response
#' is encountered, the function returns the number of the last valid page.
#'
#' Note: This approach involves multiple sequential API calls, which might be
#' slow depending on the number of pages and the rate limits of the API. If the
#' API has rate limits, pauses may be required between the calls to avoid hitting
#' those limits.
#'
#' @param base_url A character string representing the base URL of the PanelApp
#'        API endpoint. It should end with the page query parameter,
#'        e.g., "https://panelapp.genomicsengland.co.uk/api/v1/panels/?format=json&page="
#'
#' @return An integer representing the number of the last valid page.
#'
#' @examples
#' \dontrun{
#'   base_url <- "https://panelapp.genomicsengland.co.uk/api/v1/panels/?format=json&page="
#'   last_page <- find_last_page(base_url)
#'   print(last_page)
#' }
#'
#' @export
find_last_page <- function(base_url) {
  page_num <- 1
  while(TRUE) {
    response <- GET(paste0(base_url, page_num))
    content <- rawToChar(response$content)
    json_content <- fromJSON(content, flatten = TRUE)

    # Check if the response has the "detail" key indicating an invalid page
    if ("detail" %in% names(json_content) && json_content$detail == "Invalid page.") {
      return(page_num - 1)
    }
    page_num <- page_num + 1
  }
}