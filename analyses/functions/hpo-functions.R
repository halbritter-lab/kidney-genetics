require(jsonlite)
require(tidyverse)

#### This file holds analyses functions for HPO request

#' Retrieve HPO name from term ID
#'
#' This function retrieves the HPO name corresponding to a given HPO term ID.
#'
#' @param term_input_id The HPO term ID for which to retrieve the HPO name.
#'
#' @importFrom jsonlite fromJSON
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom tidyr as_tibble
#'
#' @return A tibble with the HPO name corresponding to the input HPO term ID.
#'
#' @examples
#' HPO_name_from_term("HP:1234567")
#'
#' @export
HPO_name_from_term <- function(term_input_id) {
  retries <- 3  # Number of retries before giving up
  retry_delay <- 5  # Delay in seconds between retries

  for (attempt in 1:retries) {
    tryCatch({
      hpo_term_response <- jsonlite::fromJSON(
        paste0("https://hpo.jax.org/api/hpo/term/",
               URLencode(term_input_id, reserved = TRUE)))

      hpo_term_name <- tidyr::as_tibble(hpo_term_response$details$name) %>%
        dplyr::select(hpo_mode_of_inheritance_term_name = value)

      return(hpo_term_name)
    }, error = function(e) {
      if (attempt < retries) {
        message("Retrying after error:", conditionMessage(e))
        Sys.sleep(retry_delay)
      } else {
        stop("Failed after multiple retries:", conditionMessage(e))
      }
    })
  }
}


#' Retrieve HPO definition from term ID
#'
#' This function retrieves the HPO definition corresponding to a given HPO term ID.
#'
#' @param term_input_id The HPO term ID for which to retrieve the HPO definition.
#'
#' @importFrom jsonlite fromJSON
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom tidyr as_tibble
#'
#' @return A tibble with the HPO definition corresponding to the input HPO term ID.
#'
#' @examples
#' HPO_definition_from_term("HP:1234567")
#'
#' @export
HPO_definition_from_term <- function(term_input_id) {
  retries <- 3  # Number of retries before giving up
  retry_delay <- 5  # Delay in seconds between retries

  for (attempt in 1:retries) {
    tryCatch({
      hpo_term_response <- jsonlite::fromJSON(
        paste0("https://hpo.jax.org/api/hpo/term/",
               URLencode(term_input_id, reserved = TRUE)))

      hpo_term_definition <- tidyr::as_tibble(hpo_term_response$details$definition) %>%
        dplyr::select(hpo_mode_of_inheritance_term_definition = value)

      return(hpo_term_definition)
    }, error = function(e) {
      if (attempt < retries) {
        message("Retrying after error:", conditionMessage(e))
        Sys.sleep(retry_delay)
      } else {
        stop("Failed after multiple retries:", conditionMessage(e))
      }
    })
  }
}


#' Retrieve count of HPO children from term ID
#'
#' This function retrieves the count of HPO children terms for a given HPO term ID.
#'
#' @param term_input_id The HPO term ID for which to retrieve the count of children.
#'
#' @importFrom jsonlite fromJSON
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom tidyr as_tibble
#'
#' @return An integer representing the count of HPO children terms.
#'
#' @examples
#' HPO_children_count_from_term("HP:1234567")
#'
#' @export
HPO_children_count_from_term <- function(term_input_id) {
  retries <- 3  # Number of retries before giving up
  retry_delay <- 5  # Delay in seconds between retries

  for (attempt in 1:retries) {
    tryCatch({
      hpo_term_response <- jsonlite::fromJSON(
        paste0("https://hpo.jax.org/api/hpo/term/",
               URLencode(term_input_id, reserved = TRUE)))

      hpo_term_children_count <- length(hpo_term_response$relations$children)

      return(hpo_term_children_count)
    }, error = function(e) {
      if (attempt < retries) {
        message("Retrying after error:", conditionMessage(e))
        Sys.sleep(retry_delay)
      } else {
        stop("Failed after multiple retries:", conditionMessage(e))
      }
    })
  }
}


#' This function retrieves the HPO children terms for a given HPO term ID.
#'
#' @param term_input_id The HPO term ID for which to retrieve the children terms.
#'
#' @importFrom jsonlite fromJSON
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom tidyr as_tibble
#'
#' @return A tibble with the HPO children terms corresponding to the input HPO term ID.
#'
#' @examples
#' HPO_children_from_term("HP:1234567")
#'
#' @export
HPO_children_from_term <- function(term_input_id) {
  retries <- 3  # Number of retries before giving up
  retry_delay <- 5  # Delay in seconds between retries

  for (attempt in 1:retries) {
    tryCatch({
      hpo_term_response <- jsonlite::fromJSON(
        paste0("https://hpo.jax.org/api/hpo/term/",
               URLencode(term_input_id, reserved = TRUE)))

      hpo_term_children <- tidyr::as_tibble(hpo_term_response$relations$children)

      return(hpo_term_children)
    }, error = function(e) {
      if (attempt < retries) {
        message("Retrying after error:", conditionMessage(e))
        Sys.sleep(retry_delay)
      } else {
        stop("Failed after multiple retries:", conditionMessage(e))
      }
    })
  }
}



#' Retrieve all HPO descendants and the term itself from term ID
#'
#' This function retrieves all the HPO children terms, including nested children,
#' for a given HPO term ID.
#'
#' @param term_input The HPO term ID for which to retrieve all children terms.
#' @param all_children_list A list to accumulate all children terms. Should be empty at the initial call.
#'
#' @return A tibble with all the HPO children terms, including nested children,
#'   corresponding to the input HPO term ID.
#'
#' @examples
#' HPO_all_children_from_term("HP:1234567", list())
#'
#' @export
HPO_all_children_from_term <- function(term_input, all_children_list = list()) {

  children_list <- HPO_children_from_term(term_input)

  # Combine all_children_list and term_input
  all_children_list <- c(all_children_list, term_input)

  if (length(children_list) != 0) {
    for (p in children_list$ontologyId) {
        # Update all_children_list with each recursive call
        all_children_list <- HPO_all_children_from_term(p, all_children_list)
    }
  }

  # If no more children are found, convert the list to a tibble and remove duplicates
  if (length(children_list) == 0) {
    all_children_tibble <- tidyr::as_tibble(unlist(all_children_list)) %>%
      unique()
    return(all_children_tibble)
  } else {
    return(all_children_list)
  }
}