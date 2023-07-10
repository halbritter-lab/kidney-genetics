require(tidyverse)
require(readr)
require(stringr)
require(fs)
require(lubridate)


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


#' Check the age of the most recent file in a directory
#'
#' This function checks the age of the most recent file with a given basename in a
#' specified directory. It returns TRUE if the newest file is younger than the
#' specified duration (in months), and FALSE otherwise.
#'
#' @param file_basename A string. The basename of the files to check.
#' This should be in the format "filename.", e.g. "hpo_list_kidney."
#'
#' @param folder A string. The directory where the files are located.
#'
#' @param months A numeric. The number of months to compare the file's age with.
#'
#' @return A logical. Returns TRUE if the most recent file is younger than the
#' specified number of months, and FALSE otherwise.
#'
#' @examples
#' \dontrun{
#' check_file_age("hpo_list_kidney.", "shared/", 1)
#' }
#'
#' @importFrom fs dir_ls
#' @importFrom stringr str_extract
#' @importFrom lubridate as.Date interval months
#'
#' @export
check_file_age <- function(file_basename, folder, months) {

  # Construct the regex pattern for the files
  pattern <- paste0(file_basename, "\\d{4}-\\d{2}-\\d{2}\\.csv\\.gz$")

  # Get the list of files
  files <- dir_ls(folder, regexp = pattern)

  # If there are no files, we set the time to the start of Unix epoch 
  if (length(files) == 0) {
    newest_date <- as.Date("1970-01-01")
  } else {
    # Extract the dates from the file names
    dates <- str_extract(files, "\\d{4}-\\d{2}-\\d{2}")

    # Convert the dates to Date objects
    dates <- as.Date(dates)

    # Get the newest date
    newest_date <- max(dates, na.rm = TRUE)
  }

  # Get the current date
  current_date <- Sys.Date()

  # Compute the difference in months between the current time and the newest file time
  time_diff <- interval(newest_date, current_date) / months(1)

  # Return TRUE if the newest file is older than the specified number of months, and FALSE otherwise
  return(time_diff < months)
}