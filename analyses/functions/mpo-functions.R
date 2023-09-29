require(ontologyIndex)


#' Retrieve the MP children terms for a given MP term ID.
#'
#' This function retrieves the children terms of a given term from the Mammalian
#' Phenotype Ontology (MPO) using the ontologyIndex package.
#'
#' @param term_input_id The MP term ID for which to retrieve the children terms.
#'
#' @importFrom ontologyIndex get_term_property
#' @importFrom tidyr as_tibble
#'
#' @return A tibble with the MP children terms corresponding to the input MP term ID.
#'
#' @examples
#' \dontrun{
#' # Load the Mammalian Phenotype Ontology (MPO) before using the function
#' # Example usage of the function:
#' mp_children_from_term("MP:0005502")
#' }
#'
#' @export
mp_children_from_term <- function(term_input_id) {
  # Use get_term_property to retrieve the children of the given term
  children_terms <- get_term_property(ontology = mpo, property = "children", term = term_input_id)

  # Convert the list of children terms to a tibble
  children_tibble <- tidyr::as_tibble(data.frame(ontologyId = children_terms))

  return(children_tibble)
}


#' Retrieve all MP descendants and the term itself from term ID.
#'
#' This function retrieves all the MP children terms, including nested children,
#' for a given MP term ID using the Mammalian Phenotype Ontology (MPO).
#'
#' @param term_input_id The MP term ID for which to retrieve all children terms.
#' @param all_children_list A list to accumulate all children terms. Should be empty at the initial call.
#'
#' @importFrom tidyr as_tibble
#'
#' @return A tibble with all the MP children terms, including nested children,
#'   corresponding to the input MP term ID.
#'
#' @examples
#' \dontrun{
#' # Load the Mammalian Phenotype Ontology (MPO) before using the function
#' # Example usage of the function with an empty list as the initial call:
#' mp_all_children_from_term("MP:0005502", list())
#' }
#'
#' @export
mp_all_children_from_term <- function(term_input_id, all_children_list = list()) {

  # Retrieve the immediate children of the current term
  children_tibble <- mp_children_from_term(term_input_id)

  # Combine all_children_list and term_input_id
  all_children_list <- c(all_children_list, list(term_input_id))

  if (nrow(children_tibble) != 0) {
    for (child_id in children_tibble$ontologyId) {
      # Update all_children_list with each recursive call
      all_children_list <- mp_all_children_from_term(child_id, all_children_list)
    }
  }

  # If no more children are found, flatten the list and convert it to a tibble and remove duplicates
  if (nrow(children_tibble) == 0) {
    all_children_vector <- unlist(all_children_list, use.names = FALSE)
    all_children_tibble <- tidyr::as_tibble(data.frame(ontologyId = unique(all_children_vector)))
    return(all_children_tibble)
  } else {
    return(all_children_list)
  }
}
