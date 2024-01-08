
# function that returns a list of sublists containing the subclusters of the STRING db cluster analysis
# if a subcluster is larger than 'min_number' it will divided into more subclusters
get_STRING_clusters <- function(STRING_id_vec, min_number, algorithm = "walktrap") {
  result_list <- string_db_full$get_clusters(STRING_id_vec,
                                             algorithm = algorithm)
  subgroups <- list()
  
  for (i in 1:length(result_list)) {
    if (length(result_list[[i]]) > min_number) {
      subgroups[[i]] <- get_STRING_clusters(result_list[[i]], min_number)
    } else {
      subgroups[[i]] <- result_list[[i]]
    }
  }
  return(subgroups)
}


# function that returns the index of the list in which the string of interest occurs 
get_sublist_index <- function(cluster_list, searchterm){
  if (!(searchterm %in% unlist(cluster_list))){
    return(NA)
  }
  else {
    for (i in 1:length(cluster_list)){
      if (searchterm %in% unlist(cluster_list[i])){
        return(i)
      }
    }
  }
}


# function that returns the full index of the string of interest within the list of sublists
get_full_index <- function(cluster_list, searchterm, res = c()){
  if (!(searchterm %in% unlist(cluster_list))){
    return(NA)
  }
  else if (length(cluster_list) == 1){
    index <- which(unlist(cluster_list) == searchterm)
    full_index <- c(res, index)
    return(paste(full_index[-length(full_index)], collapse = "-"))
  }
  else{
    sublist_index <- get_sublist_index(cluster_list, searchterm=searchterm)
    res <- c(res, sublist_index)
    get_full_index(cluster_list[[sublist_index]], searchterm, res=res)
  }
} 


# function to plot a pie chart showing the kidney disease group distribution within the subcluster
plot_disease_group_distribution <- function(subcluster, disease_group_df){
  
  
  # # define a custom color palette for kidney_disease_group_short
  # custom_colors <- c("tubulopathy" = "Red", 
  #                    "glomerulopathy" = "Green", 
  #                    "cancer" = "Blue", 
  #                    "cakut" = "Purple", 
  #                    "cyst_cilio" = "Orange", 
  #                    "complement" = "Yellow",
  #                    "nephrocalcinosis" = "Grey") 
  
  # define a custom (colorblind readable) color palette for kidney_disease_group_short
  custom_colors <- c("tubulopathy" = "#88CCEE", 
                     "glomerulopathy" = "#CC6677", 
                     "cancer" = "#DDCC77", 
                     "cakut" = "#44AA99", 
                     "cyst_cilio" = "#AA4499", 
                     "complement" = "#999933",
                     "nephrocalcinosis" = "#888888")
  
  
  # get subset df based on specified subcluster
  subset_df <- disease_group_df %>% 
    filter(grepl(paste0("^", subcluster), cluster_index))
  
  if (nrow(subset_df) == 0){
    message("Subcluster ", subcluster, " does not exist.")
  }
  
  # plot the pie chart
  pie_chart <- ggplot(subset_df, aes(x = "", fill = kidney_disease_group_short)) +
    geom_bar(width = 1) +
    scale_fill_manual(values = custom_colors) +
    coord_polar(theta = "y") +
    labs(title = paste("Full network - Cluster",  subcluster, "| Total Instances:", nrow(subset_df)))
  
  # define the filename
  # filename <- paste0("results/prot_interact_full_network_cluster_", subcluster, "_pie_chart.", current_date, ".png")
  
  # save the pie chart as a PNG file
  # ggsave(filename, plot = pie_chart, width = 5, height = 5) # TODO: where to store plots?
  
  return(pie_chart)
}


# function to get all interactions above a specified mininum STRING combined_score
get_all_interactions_above_score <- function(STRING_id_vec, string_db, min_comb_score){
  # get all interactions above minimum combined score
  all_interactions <- string_db_full$get_interactions(STRING_id_vec) %>% 
    distinct() %>% 
    filter(combined_score >= min_comb_score)
  
  # modify df so that each gene pair appears only once (TODO: can possibly be omitted)
  all_interactions <- all_interactions %>%  
    rowwise %>%
    mutate(new_from = sort(c(from, to))[1],
           new_to = sort(c(from, to))[2]) %>% 
    dplyr::select(-from, -to) %>% 
    rename(from = new_from, to = new_to) %>% distinct()
  
  # create a column that contains both involved genes per row, respectively
  all_interactions <- all_interactions %>% 
    mutate(both_genes = list(c(from, to))) %>% 
    mutate(row_index = cur_group_id())
  
  return(all_interactions)
}



# function to get direct contacts of a gene and remove corresponding rows in connection df
get_direct_contacts <- function(index_gene, connection_df, min_comb_score){
  # get all connections of the index gene above a minimum combined score
  gene_connections <- connection_df %>% 
    filter(index_gene %in% both_genes, combined_score >= min_comb_score)
  
  # remove index gene itself from direct contacts vector
  direct_contacts <- setdiff(union(gene_connections$from, gene_connections$to), index_gene)
  
  # remove rows from connection df that have already been checked here (for further analysis)
  sub_connection_df <- connection_df %>% filter(!(row_index %in% gene_connections$row_index))
  
  return(list(direct_contacts = direct_contacts, sub_connection_df = sub_connection_df))
}


# function to all direct and indirect contacts of an index gene (without recursive call, faster)
get_all_contacts <- function(index_gene, connection_df, min_comb_score){
  # get direct contacts of index_gene
  con_list <- get_direct_contacts(index_gene, connection_df, min_comb_score)
  
  # create a results vector
  res <- c(index_gene)
  
  # create a vector for genes that have to be checked for further contacts
  genes_to_check <- con_list$direct_contacts
  
  # check all genes in genes_to_check for further contacts
  while (length(genes_to_check) > 0){
    gene1 <- genes_to_check[1]
    
    # get direct contacts of gene that is checked
    con_list <- get_direct_contacts(gene1, con_list$sub_connection_df, min_comb_score)
    
    # add gene that is checked to results vector
    res <- union(res, gene1)
    
    # removed checked gene from genes_to_check
    genes_to_check <- setdiff(genes_to_check, gene1)
    
    # add new genes (direct contacts of checked gene) to genes_to_check
    genes_to_check <- union(genes_to_check, setdiff(con_list$direct_contacts, res))
  }
  
  return(res)
}




# function to create an edglist from all direct and indirect contacts of an index gene
# (and the index gene itself) above a minimum combined STRING score
create_edgelist <- function(connection_df, all_contacts_to_index, min_comb_score, symbol_annotation_df){
  # filter connection df for genes to plot
  interactions <- connection_df %>% 
    filter(from %in% all_contacts_to_index, to %in% all_contacts_to_index, combined_score > min_comb_score)
  
  # annotate with symbols
  interactions <- interactions %>% 
    left_join(symbol_annotation_df, by = c("from" = "STRING_id")) %>% 
    rename(symbol_from = symbol) %>% 
    left_join(symbol_annotation_df, by = c("to" = "STRING_id")) %>% 
    rename(symbol_to = symbol) %>% 
    distinct() %>% 
    dplyr::select(symbol_from, symbol_to, combined_score)  
  
  # create an edgelist
  edgel <- cbind(interactions$symbol_from, interactions$symbol_to)
  
  return(edgel)
}


# function to plot the interaction network given by an edgelist
plot_interaction_network <- function(edgelist, disease_group_df){
  
  # create a network from the edgelist
  netw <- network(edgelist, directed = F)
  
  # get unique vertex names
  unique_vertices <- unique(c(edgelist[, 1], edgelist[, 2]))
  
  # define a custom (colorblind readable) color palette for kidney_disease_group_short
  custom_colors <- c("tubulopathy" = "#88CCEE", 
                     "glomerulopathy" = "#CC6677", 
                     "cancer" = "#DDCC77", 
                     "cakut" = "#44AA99", 
                     "cyst_cilio" = "#AA4499", 
                     "complement" = "#999933",
                     "nephrocalcinosis" = "#888888")
  
  # define vertex colors
  vertex_colors <- disease_group_df %>% 
    filter(symbol %in% unique_vertices) %>% 
    mutate(color = custom_colors[kidney_disease_group_short])
  
  vertex_color_map <- setNames(vertex_colors$color, vertex_colors$symbol)
  
  # assign colors to vertices based on their names
  vertex_colors <- sapply(network.vertex.names(netw), function(v) vertex_color_map[[v]])
  
  # create a plotly object
  p <- ggnet2(netw, label = TRUE, color = vertex_colors, label.size = 2.5)
  
  p_plotly <- ggplotly(p)
  
  # add legend to the plot
  p_plotly <- p_plotly %>%
    add_trace(
      x = rep(-0.1, length(custom_colors)),
      y = seq(-0.2, length.out = length(custom_colors), by = -0.04),
      mode = "markers+text",
      marker = list(color = unname(custom_colors)), 
      text = names(custom_colors),
      marker = list(size = 10),
      showlegend = TRUE, 
      legendgroup = "vertex_colors",
      textposition = "middle right",
      textfont = list(size = 11, color = "black")
    ) 
  
  # remove axis ticks and tick labels
  p_plotly <- p_plotly %>%
    layout(
      xaxis = list(
        showticklabels = FALSE,
        showline = FALSE,
        showgrid = FALSE,
        zeroline = FALSE
      ),
      yaxis = list(
        showticklabels = FALSE,
        showline = FALSE,
        showgrid = FALSE,
        zeroline = FALSE
      )
    )
  
  # disable hover info
  p_plotly <- style(p_plotly, hoverinfo = "none") 
  
  return(p_plotly)
}

# function to create a plotly plot of a network of direct and indirect contacts with a of an index gene
# minimum combined score in STRING
plot_network_of_index_gene <- function(index_gene, string_db, min_comb_score, STRING_id_vec, disease_group_df){
  
  # get all STRING interactions between all genes in STRING_id_vec with a minimum combined score in STRING
  all_interactions <- get_all_interactions_above_score(STRING_id_vec = STRING_id_vec,
                                                       string_db = string_db,
                                                       min_comb_score = min_comb_score)
  # get all contacts of the index_gene above the minimum combined score in STRING
  
  all_contacts_to_index <- get_all_contacts(index_gene = index_gene,
                                            connection_df = all_interactions,
                                            min_comb_score = min_comb_score)
  
  # create an edgelist of all contacts of the index gene (direct/indirect contacts)
  edgelist <- create_edgelist(connection_df = all_interactions, 
                              all_contacts_to_index = all_contacts_to_index, 
                              min_comb_score = min_comb_score, 
                              symbol_annotation_df = distinct(disease_group_df[c("STRING_id", "symbol")]))
  
  # plot interaction network of index gene
  interaction_plot <- plot_interaction_network(edgelist = edgelist,
                                               disease_group_df = disease_group_df)
  
  return(interaction_plot)
}



