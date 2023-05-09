#### This file holds analyses functions


## functions for HPO request

HPO_name_from_term <- function(term_input_id) {
	hpo_term_response <- fromJSON(paste0("https://hpo.jax.org/api/hpo/term/", URLencode(term_input_id, reserved=T)))
	hpo_term_name <- as_tibble(hpo_term_response$details$name) %>%
	select(hpo_mode_of_inheritance_term_name = value)

	return(hpo_term_name)
}


HPO_definition_from_term <- function(term_input_id) {
	hpo_term_response <- fromJSON(paste0("https://hpo.jax.org/api/hpo/term/", URLencode(term_input_id, reserved=T)))
	hpo_term_definition <- as_tibble(hpo_term_response$details$definition) %>%
	select(hpo_mode_of_inheritance_term_definition = value)

	return(hpo_term_definition)
}


HPO_children_count_from_term <- function(term_input_id) {
	hpo_term_response <- fromJSON(paste0("https://hpo.jax.org/api/hpo/term/", URLencode(term_input_id, reserved=T)))
	hpo_term_children_count <- as_tibble(hpo_term_response$relations$children)

	return(length(hpo_term_children_count))
}


HPO_children_from_term <- function(term_input_id) {
	hpo_term_response <- fromJSON(paste0("https://hpo.jax.org/api/hpo/term/", URLencode(term_input_id, reserved=T)))
	hpo_term_children <- as_tibble(hpo_term_response$relations$children)

	return(hpo_term_children)
}

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