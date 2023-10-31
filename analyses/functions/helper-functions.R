

require(dplyr)
require(tibble)

#' Percentile Normalization of a Column in a Tibble
#'
#' This function normalizes a specified column in a tibble by transforming each value
#' to the percentile in which it falls relative to the entire column. 
#'
#' @param data A tibble containing the data.
#' @param colname A string representing the name of the column to be normalized.
#'
#' @return A tibble with an additional column named "<colname>_percentile" where each
#'         entry represents the percentile of the corresponding value in the original
#'         column. Percentile values range from 0 to 1.
#'
#' @examples
#' library(tibble)
#' library(dplyr)
#'
#' # Example tibble
#' df <- tibble::tibble(
#'   approved_symbol = c("ABCG2", "ACE", "ACTA2"),
#'   source_count = c(3, 10, 6)
#' )
#'
#' # Applying normalization
#' df_normalized <- normalize_percentile(df, "source_count")
#'
#' @export
normalize_percentile <- function(data, colname) {
  data %>%
    mutate("{colname}_percentile" := rank(!!sym(colname), ties.method = "average") / n()) 
}


#' Get Current Date in ISO 8601 Format
#'
#' This function returns the current date in ISO 8601 format ("YYYY-MM-DD"). 
#' The date is calculated based on Coordinated Universal Time (UTC).
#'
#' @return A character string representing the current date in "YYYY-MM-DD" format.
#'
#' @examples
#' current_date <- get_current_date_iso8601()
#' print(current_date)
#'
#' @export
get_current_date_iso8601 <- function() {
  strftime(as.POSIXlt(Sys.time(), "UTC"), "%Y-%m-%d")
}