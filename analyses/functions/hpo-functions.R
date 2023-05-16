#### This file holds analyses functions for HPO request

#' Retrieve HPO name from term ID
#'
#' This function retrieves the HPO name corresponding to a given HPO term ID.
#'
#' @param term_input_id The HPO term ID for which to retrieve the HPO name.
#'
#' @return A tibble with the HPO name corresponding to the input HPO term ID.
#'
#' @examples
#' HPO_name_from_term("HPO:1234567")
#'
#' @export
HPO_name_from_term <- function(term_input_id) {
  hpo_term_response <- fromJSON(paste0("https://hpo.jax.org/api/hpo/term/", URLencode(term_input_id, reserved=T)))
  hpo_term_name <- as_tibble(hpo_term_response$details$name) %>%
  select(hpo_mode_of_inheritance_term_name = value)

  return(hpo_term_name)
}


#' Retrieve HPO definition from term ID
#'
#' This function retrieves the HPO definition corresponding to a given HPO term ID.
#'
#' @param term_input_id The HPO term ID for which to retrieve the HPO definition.
#'
#' @return A tibble with the HPO definition corresponding to the input HPO term ID.
#'
#' @examples
#' HPO_definition_from_term("HPO:1234567")
#'
#' @export
HPO_definition_from_term <- function(term_input_id) {
  hpo_term_response <- fromJSON(paste0("https://hpo.jax.org/api/hpo/term/", URLencode(term_input_id, reserved=T)))
  hpo_term_definition <- as_tibble(hpo_term_response$details$definition) %>%
  select(hpo_mode_of_inheritance_term_definition = value)

  return(hpo_term_definition)
}


#' Retrieve count of HPO children from term ID
#'
#' This function retrieves the count of HPO children terms for a given HPO term ID.
#'
#' @param term_input_id The HPO term ID for which to retrieve the count of children.
#'
#' @return An integer representing the count of HPO children terms.
#'
#' @examples
#' HPO_children_count_from_term("HPO:1234567")
#'
#' @export
HPO_children_count_from_term <- function(term_input_id) {
  hpo_term_response <- fromJSON(paste0("https://hpo.jax.org/api/hpo/term/", URLencode(term_input_id, reserved=T)))
  hpo_term_children_count <- as_tibble(hpo_term_response$relations$children)

  return(length(hpo_term_children_count))
}


#' This function retrieves the HPO children terms for a given HPO term ID.
#'
#' @param term_input_id The HPO term ID for which to retrieve the children terms.
#'
#' @return A tibble with the HPO children terms corresponding to the input HPO term ID.
#'
#' @examples
#' HPO_children_from_term("HPO:1234567")
#'
#' @export
HPO_children_from_term <- function(term_input_id) {
  hpo_term_response <- fromJSON(paste0("https://hpo.jax.org/api/hpo/term/", URLencode(term_input_id, reserved=T)))
  hpo_term_children <- as_tibble(hpo_term_response$relations$children)

  return(hpo_term_children)
}


#' Retrieve all HPO children from term ID
#'
#' This function retrieves all the HPO children terms, including nested children,
#' for a given HPO term ID.
#'
#' @param term_input The HPO term ID for which to retrieve all children terms.
#'
#' @return A tibble with all the HPO children terms, including nested children,
#'   corresponding to the input HPO term ID.
#'
#' @examples
#' HPO_all_children_from_term("HPO:1234567")
#'
#' @export
HPO_all_children_from_term <- function(term_input) {

  children_list <- HPO_children_from_term(term_input)
  all_children_list <<- append(all_children_list, term_input)

  if(length(children_list)!=0)
  {
    for (p in children_list$ontologyId) {
        all_children_list <<- append(all_children_list, p)
        Recall(p)
    }
  }
  all_children_tibble <- as_tibble(unlist(all_children_list)) %>% unique

  return(all_children_tibble)
}