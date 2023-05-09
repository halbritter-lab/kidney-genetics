#### This file holds analyses functions


## functions for NCBI datasets v2 API request

# this function uses the NCBI datasets v2 API to return gene information based on a NCBI gene ID
gene_info_from_gene_id <- function(input, api_key = NCBI_API_KEY, request_max = 20) {
	ep_base_url <- "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/gene/id/"
	ep_fields_args <- "?table_fields=gene-id&table_fields=gene-type&table_fields=description"

	# TODO: implement retry if empty	
	# TODO: implement error handling if input not a number or NA

	input_list <- as_tibble(input, .name_repair = "unique") %>%
		mutate(value = as.character(value))
	input_list_unique <- as_tibble(input, .name_repair = "unique") %>%
		unique()

	row_number <- nrow(input_list_unique)
	groups_number <- ceiling(row_number/request_max)
	
	input_list_result <- input_list_unique %>%
		mutate(group = sample(1:groups_number, row_number, replace=T)) %>%
		group_by(group) %>%
		summarise(query = paste(unique(value), collapse = ",")) %>%
		ungroup() %>%
		mutate(url = paste0(ep_base_url,
			query,
			ep_fields_args)) %>%
		rowwise() %>%
		mutate(results = list(fromJSON(content(GET(url, add_headers("api-key" = api_key), accept_json()), "text", encoding = "UTF-8"))$reports$gene %>% 
			as_tibble() %>%
			unnest(annotations) %>%
			unnest(nomenclature_authority) %>%
			unnest(ensembl_gene_ids, keep_empty = TRUE) %>%
			select(gene_id, symbol, tax_id, taxname, , identifier, ensembl_gene_ids) %>%
			unique())) %>%
		select(-group, - query, -url) %>%
		unnest(results)

	input_list_return <- input_list %>%
		left_join(input_list_result, by = c("value" = "gene_id"))

	return(input_list_return)
}

# this function finds genes from a pubtator API request
pubtator_genes_in_request <- function(query, page, filter_type = "Gene") {
	url <- paste0("https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/search?q=", query, "&page=", page)
	search_request <- fromJSON(URLencode(url), flatten = TRUE)

	# TODO: implement retry if error
	# TODO: also return NCBI gene ID

	search_results_tibble <- search_request$results %>%
		as_tibble()
		
	search_results_filtered <- search_results_tibble %>%
		{if( !("accessions" %in% colnames(.)) ) add_column(., accessions = "NULL") else .} %>%
		filter(accessions != "NULL")

	if (nrow(search_results_filtered) > 0)
	{
		search_results <- search_results_filtered %>%
			select(pmid, passages) %>%
			unnest(passages) %>%
			select(pmid, text_part = infons.type, annotations) %>%
			rowwise() %>%
			mutate(empty = is_empty(annotations)) %>%
			ungroup() %>%
			filter(annotations != "NULL" & !empty) %>%
			unnest(annotations, names_repair = "universal") %>%
			select(pmid, text_part, text, type = infons.type, text_identifier = infons.identifier) %>%
			unique() %>%
			filter(type == filter_type) %>%
			group_by(pmid, text, text_identifier) %>%
			summarise(text_part = paste(unique(text_part), collapse = " | "),
				source = paste(unique(type), collapse = " | "),
				text_part_count = n(),
				.groups = "keep") %>%
			ungroup()

		return(search_results)
	} else {
		return(NULL)
	}
	}


# a function that returns the number of pages for a request to the NCBI PubTator API
pubtator_pages_request <- function(query) {
	url <- paste0("https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/search?q=", query, "&page=", 1)
	search_request <- fromJSON(URLencode(url), flatten = TRUE)

	return(search_request$total_pages)
}