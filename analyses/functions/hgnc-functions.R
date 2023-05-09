#### This file holds analyses functions


## functions for HGNC IDs
hgnc_id_from_prevsymbol <- function(symbol_input)  {
	symbol_request <- fromJSON(paste0("http://rest.genenames.org/search/prev_symbol/", symbol_input))

	hgnc_id_from_symbol <- as_tibble(symbol_request$response$docs)
	
	hgnc_id_from_symbol <- hgnc_id_from_symbol %>%
	mutate(hgnc_id = if (exists('hgnc_id', where = hgnc_id_from_symbol)) hgnc_id else NA) %>%
	mutate(symbol = if (exists('symbol', where = hgnc_id_from_symbol)) symbol else "") %>%
	mutate(score = if (exists('score', where = hgnc_id_from_symbol)) score else 0) %>%
	arrange(desc(score)) %>%
	mutate(hgnc_id = as.integer(str_split_fixed(hgnc_id, ":", 2)[, 2]))

	return(as.integer(hgnc_id_from_symbol$hgnc_id[1]))
}

hgnc_id_from_aliassymbol <- function(symbol_input)  {
	symbol_request <- fromJSON(paste0("http://rest.genenames.org/search/alias_symbol/", symbol_input))

	hgnc_id_from_symbol <- as_tibble(symbol_request$response$docs)
	
	hgnc_id_from_symbol <- hgnc_id_from_symbol %>%
	mutate(hgnc_id = if (exists('hgnc_id', where = hgnc_id_from_symbol)) hgnc_id else NA) %>%
	mutate(symbol = if (exists('symbol', where = hgnc_id_from_symbol)) symbol else "") %>%
	mutate(score = if (exists('score', where = hgnc_id_from_symbol)) score else 0) %>%
	arrange(desc(score)) %>%
	mutate(hgnc_id = as.integer(str_split_fixed(hgnc_id, ":", 2)[, 2]))

	return(as.integer(hgnc_id_from_symbol$hgnc_id[1]))
}

hgnc_id_from_symbol <- function(symbol_tibble) {
	symbol_list_tibble <- as_tibble(symbol_tibble) %>% select(symbol = value) %>% mutate(symbol = toupper(symbol))
	
	symbol_request <- fromJSON(paste0("http://rest.genenames.org/search/symbol/", str_c(symbol_list_tibble$symbol, collapse = "+OR+")))

	hgnc_id_from_symbol <- as_tibble(symbol_request$response$docs)
	
	hgnc_id_from_symbol <- hgnc_id_from_symbol %>%
	mutate(hgnc_id = if (exists('hgnc_id', where = hgnc_id_from_symbol)) hgnc_id else NA) %>%
	mutate(symbol = if (exists('symbol', where = hgnc_id_from_symbol)) toupper(symbol) else "") %>%
	mutate(hgnc_id = as.integer(str_split_fixed(hgnc_id, ":", 2)[, 2]))
		
	return_tibble <- symbol_list_tibble %>% 
	left_join(hgnc_id_from_symbol, by = "symbol") %>%
	select(hgnc_id)

	return(return_tibble)
}	

hgnc_id_from_symbol_grouped <- function(input_tibble, request_max = 150) {
	input_tibble <- as_tibble(input_tibble)
	
	row_number <- nrow(input_tibble)
	groups_number <- ceiling(row_number/request_max)
	
	input_tibble_request <- input_tibble %>%
	mutate(group = sample(1:groups_number, row_number, replace=T)) %>%
	group_by(group) %>%
	mutate(response = hgnc_id_from_symbol(value)$hgnc_id) %>%
	ungroup()
	
	input_tibble_request_repair <- input_tibble_request %>%
	filter(is.na(response)) %>%
	select(value) %>%
	unique() %>%
	rowwise() %>%
	mutate(response = hgnc_id_from_prevsymbol(value)) %>%
	mutate(response = case_when(!is.na(response) ~ response, is.na(response) ~ hgnc_id_from_aliassymbol(value)))
	
	input_tibble_request <- input_tibble_request %>%
	left_join(input_tibble_request_repair, by = "value") %>%
	mutate(response = case_when(!is.na(response.x) ~ response.x, is.na(response.x) ~ response.y))
	
	return(input_tibble_request$response)
}

symbol_from_hgnc_id <- function(hgnc_id_tibble) {
	hgnc_id_list_tibble <- as_tibble(hgnc_id_tibble) %>%
		select(hgnc_id = value) %>%
		mutate(hgnc_id = as.integer(hgnc_id))
	
	hgnc_id_request <- fromJSON(paste0("http://rest.genenames.org/search/hgnc_id/", str_c(hgnc_id_list_tibble$hgnc_id, collapse = "+OR+")))

	hgnc_id_from_hgnc_id <- as_tibble(hgnc_id_request$response$docs)
	
	hgnc_id_from_hgnc_id <- hgnc_id_from_hgnc_id %>%
	mutate(hgnc_id = if (exists('hgnc_id', where = hgnc_id_from_hgnc_id)) hgnc_id else NA) %>%
	mutate(hgnc_id = if (exists('hgnc_id', where = hgnc_id_from_hgnc_id)) toupper(hgnc_id) else "") %>%
	mutate(hgnc_id = as.integer(str_split_fixed(hgnc_id, ":", 2)[, 2]))
		
	return_tibble <- hgnc_id_list_tibble %>% 
	left_join(hgnc_id_from_hgnc_id, by = "hgnc_id") %>%
	select(symbol)

	return(return_tibble)
}


symbol_from_hgnc_id_grouped <- function(input_tibble, request_max = 150) {
	input_tibble <- as_tibble(input_tibble)
	
	row_number <- nrow(input_tibble)
	groups_number <- ceiling(row_number/request_max)
	
	input_tibble_request <- input_tibble %>%
		mutate(group = sample(1:groups_number, row_number, replace=T)) %>%
		group_by(group) %>%
		mutate(response = symbol_from_hgnc_id(value)$symbol) %>%
		ungroup()
	
	return(input_tibble_request$response)
}