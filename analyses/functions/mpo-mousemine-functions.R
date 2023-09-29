require(ontologyIndex)
require(InterMineR)


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
  children_terms <- get_term_property(ontology = mpo, property = "children",
                                      term = term_input_id)

  # Convert the list of children terms to a tibble
  children_tibble <- tidyr::as_tibble(data.frame(ontology_id = children_terms))

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
    for (child_id in children_tibble$ontology_id) {
      # Update all_children_list with each recursive call
      all_children_list <- mp_all_children_from_term(child_id, all_children_list)
    }
  }

  # If no more children are found, flatten the list, convert it to a tibble and remove duplicates
  if (nrow(children_tibble) == 0) {
    all_children_vector <- unlist(all_children_list, use.names = FALSE)
    all_children_tibble <- tidyr::as_tibble(data.frame(ontology_id = unique(all_children_vector)))
    return(all_children_tibble)
  } else {
    return(all_children_list)
  }
}


#' Compare Gene Phenotypes
#'
#' This function retrieves the phenotypes associated with a given gene from the
#' MouseMine database using the InterMineR package. It then compares these
#' phenotypes with a provided list of phenotype identifiers and returns the
#' comparison results grouped by zygosity.
#'
#' @param gene_name A character string specifying the name of the gene.
#' @param phenotype_ids A tibble containing a column named 'ontology_id' with
#'   the phenotype identifiers to compare against.
#'
#' @return A character string representing the comparison results, grouped by
#'   zygosity. The format of the result is "hm (true/false); ht (true/false)",
#'   where "true" indicates the presence of at least one of the terms in the
#'   input list of phenotype identifiers, and "false" indicates the absence of
#'   all terms. If there are no phenotypes or a zygosity does not exist, the
#'   result will be "NA" for that zygosity.
#'
#' @examples
#' phenotype_ids <- tibble::tibble(ontology_id = c("MP:0005502", "MP:0001756"))
#' compare_gene_phenotypes("PKD1", phenotype_ids)
#'
#' @export
compare_gene_phenotypes <- function(gene_name, phenotype_ids) {

  # Check if the required column is present in the phenotype_ids tibble
  if (!"ontology_id" %in% names(phenotype_ids)) {
    stop("The phenotype_ids tibble must contain a column named 'ontology_id'")
  }

  # Initialise the connection to the MouseMine database
  im <- initInterMine(mine = listMines()["MouseMine"])

  # Get the template query for the "GenotypePhenotype" template
  query <- getTemplateQuery(im = im, name = "_Genotype_Phenotype")

  # Set the query parameter for the gene of interest
  query$where[[4]][["value"]] <- gene_name

  # Run the query
  res <- runQuery(im, query)

  # Check if the result is empty
  if (nrow(res) == 0) {
    return("hm (NA); ht (NA)")
  }

  # Convert the result to a tibble and filter it
  res_tibble <- res %>%
    as_tibble() %>%
    filter(OntologyAnnotation.subject.zygosity %in% c("hm", "ht")) %>%
    select(
      primary_identifier = OntologyAnnotation.subject.primaryIdentifier,
      zygosity = OntologyAnnotation.subject.zygosity,
      ontology_id = OntologyAnnotation.ontologyTerm.identifier,
      ontology_name = OntologyAnnotation.ontologyTerm.name
    )

  # Compare the phenotypes and construct the result string
  hm_result <- ifelse(any(res_tibble$zygosity == "hm" & res_tibble$ontology_id %in% phenotype_ids$ontology_id), "true", "false")
  ht_result <- ifelse(any(res_tibble$zygosity == "ht" & res_tibble$ontology_id %in% phenotype_ids$ontology_id), "true", "false")

  result_string <- paste0("hm (", hm_result, "); ht (", ht_result, ")")

  return(result_string)
}
