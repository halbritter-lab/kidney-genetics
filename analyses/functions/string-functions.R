
#' Format Numeric Values as Percentage Strings
#'
#' This function takes a numeric value (or vector of values) representing a proportion
#' (e.g., 0.25 for 25%), multiplies it by 100, and then converts it into a string
#' representation with a percentage symbol.
#'
# based on: https://statisticsglobe.com/format-number-as-percentage-in-r
#'
#' @param x A numeric value or vector representing the proportion(s) to be formatted.
#' @param digits An integer specifying the number of decimal places
#'        to display in the percentage. Default is 2.
#' @param format A character string indicating the type of formatting 
#'        desired. Default is "f" for fixed. Other options include "e" for scientific notation,
#'        and "g" for general format. See \code{\link[base]{formatC}} for more details.
#' @param ... Additional arguments to be passed to the \code{\link[base]{formatC}} function.
#'
#' @return A character string or vector of strings representing the percentage(s).
#'
#' @examples
#' \dontrun{
#'   value <- 0.25
#'   percent_string <- percent(value)
#'   print(percent_string)
#'
#'   value_vector <- c(0.25, 0.50, 0.75)
#'   percent_strings <- percent(value_vector)
#'   print(percent_strings)
#' }
#'
#' @seealso \code{\link[base]{formatC}} for the underlying formatting function.
#'
#' @export
percent <- function(x, digits = 2, format = "f", ...) {
  paste0(formatC(x * 100, format = format, digits = digits, ...), "%")
}


#' Format a Number as a String with Specified Decimals and Big Mark
#'
#' This function takes a number, rounds it to the specified number of decimal places,
#' formats it with the desired big mark (e.g., comma for thousands), and then
#' converts it to a string representation.
#'
#' @param number A numeric value or a value that can be coerced into numeric.
#' @param decimals An integer specifying the number of decimal places
#'        to round to. Default is 0.
#' @param big_mark A character string specifying the symbol to use
#'        for separating every three integers. Default is ",".
#'
#' @return A character string representing the formatted number.
#'
#' @examples
#' \dontrun{
#'   value <- 1234567.89
#'   formatted_string <- format_number(value)
#'   print(formatted_string)
#'
#'   formatted_string_with_decimals <- format_number(value, decimals = 2)
#'   print(formatted_string_with_decimals)
#' }
#'
#' @export
format_number <- function(number, decimals = 0, big_mark = ",") {
  # Convert to numeric
  number <- as.numeric(number)

  # Round the number
  rounded_number <- round(number, decimals)

  # Format the rounded number
  formatted_number <- format(rounded_number, nsmall = decimals, big_mark = big_mark)

  # Convert to string
  return(toString(formatted_number))
}
