require(tidyverse) # Load the tidyverse package before using the functions
require(ontologyIndex) # Load the ontologyIndex package before using the functions
require(InterMineR) # Load the InterMineR package before using the functions

#### This file holds analyses functions for MPO request
#' Retrieve the MPO children terms for a given MPO term ID.
#'
#' This function retrieves the immediate children terms of a specified term 
#' from the Mammalian Phenotype Ontology (MPO). It uses the `ontologyIndex` package 
#' to interact with the ontology and the `tidyr` package to format the results.
#'
#' @param term_input_id The MPO term ID for which to retrieve the children terms. 
#'   This should be a valid MPO term ID string.
#'
#' @details The function starts by querying the ontology for the children of the 
#'   given term using `get_term_property()`. The result is then transformed into 
#'   a tibble. Additionally, the current system date is added as an attribute 
#'   to each term, representing the date of the query.
#'
#' @importFrom ontologyIndex get_term_property
#' @importFrom tidyr as_tibble
#'
#' @return A tibble with two columns: 
#'   - `term`: Contains the MPO children terms corresponding to the input MPO term ID.
#'   - `query_date`: Contains the date when the query was made.
#'
#' @examples
#' \dontrun{
#'   # Load the Mammalian Phenotype Ontology (MPO) before using the function
#'   # Example usage of the function:
#'   mpo_children_from_term("MP:0005502")
#' }
#'
#' @export
mpo_children_from_term <- function(term_input_id) {
  # Query the ontology for children of the given term
  children_terms <- get_term_property(ontology = mpo, property = "children", term = term_input_id)

  # Convert the result into a tibble and add the current query date
  children_tibble <- children_terms %>%
    tidyr::as_tibble() %>%
    mutate(query_date = Sys.Date()) %>%
    select(term = value, query_date)

  return(children_tibble)
}


#' Retrieve all MPO descendants and the term itself from term ID.
#'
#' @param term_input_id The MPO term ID for which to retrieve all children terms.
#' @param all_children_list A list to accumulate all children terms. Initial call should be empty.
#' @return A tibble with unique MPO children terms and the query date, including nested children,
#'   corresponding to the input MPO term ID.
#' @examples
#' \dontrun{
#' # Load the Human Phenotype Ontology (MPO) before using the function
#' # Example usage of the function with an empty list as the initial call:
#' mpo_all_children_from_term("MP:0005502")
#' }
#' @export
mpo_all_children_from_term <- function(term_input_id, all_children_list = tibble()) {

  # Retrieve the immediate children of the current term
  children_tibble <- mpo_children_from_term(term_input_id)

  # Add the current term to the results
  all_children_list <- bind_rows(all_children_list, tibble(term = term_input_id, query_date = Sys.Date()))

  # If there are children, recursively get their children
  if (nrow(children_tibble) > 0) {
    for (child_id in children_tibble$term) {
      all_children_list <- mpo_all_children_from_term(child_id, all_children_list)
    }
  }

  # Return only unique rows
  return(distinct(all_children_list))
}


#' Compare Gene Phenotypes
#'
#' This function retrieves the phenotypes associated with a given gene from the
#' MouseMine database using the InterMineR package. It then compares these
#' phenotypes with a provided list of phenotype identifiers and returns the
#' comparison results grouped by zygosity.
#'
#' @param gene_name A character string specifying the name of the gene.
#' @param phenotype_ids A tibble containing a column named 'term' with
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
#' phenotype_ids <- tibble::tibble(term = c("MP:0005502", "MP:0001756"))
#' compare_gene_phenotypes("PKD1", phenotype_ids)
#'
#' @export
compare_gene_phenotypes <- function(gene_name, phenotype_ids) {

  # Check if the required column is present in the phenotype_ids tibble
  if (!"term" %in% names(phenotype_ids)) {
    stop("The phenotype_ids tibble must contain a column named 'term'")
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
      term = OntologyAnnotation.ontologyTerm.identifier,
      ontology_name = OntologyAnnotation.ontologyTerm.name
    )

  # Compare the phenotypes and construct the result string
  hm_result <- ifelse(any(res_tibble$zygosity == "hm" & res_tibble$term %in% phenotype_ids$term), "true", "false")
  ht_result <- ifelse(any(res_tibble$zygosity == "ht" & res_tibble$term %in% phenotype_ids$term), "true", "false")

  result_string <- paste0("hm (", hm_result, "); ht (", ht_result, ")")

  return(result_string)
}
