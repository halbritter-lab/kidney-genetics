require(tidyverse)
require(readr)
require(stringr)

#' Replace multiple strings in a text file
#'
#' This function takes in the path to a text file and two vectors: one of strings
#' to find and one of strings to replace them with. It replaces each instance of
#' each string to find with the corresponding string to replace it with, and writes
#' the result to a new text file.
#'
#' @param input_file A character string specifying the path to the text file to
#'   modify.
#' @param output_file A character string specifying the path to write the modified
#'   text file to.
#' @param find_vector A character vector of strings to find in the text file.
#' @param replace_vector A character vector of strings to replace the found strings
#'   with. Must be the same length as find_vector.
#'
#' @importFrom readr read_file write_file
#' @importFrom stringr str_replace_all
#' @importFrom magrittr %>%
#'
#' @return Writes the modified text to a file. No value is returned.
#' @export
#'
#' @examples
#' \dontrun{
#' replace_strings("input.txt", "output.txt",
#'   c("find_this", "find_that"),
#'   c("replace_with_this", "replace_with_that"))
#' }
replace_strings <- function(input_file, output_file, find_vector, replace_vector) {
  # Read the file as a character vector
  text <- read_file(input_file)

  # Check that the find and replace vectors are the same length
  if (length(find_vector) != length(replace_vector)) {
    stop("Find and replace vectors must be the same length")
  }

  # Replace the strings
  new_text <- text
  for (i in seq_along(find_vector)) {
    new_text <- str_replace_all(new_text, pattern = find_vector[i],
        replacement = replace_vector[i])
  }

  # Write the new text to a file
  write_file(new_text, output_file)
}
