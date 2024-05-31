
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


# function to get all direct and indirect contacts of an index gene (without recursive call, faster)
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


# function to get all direct and indirect contacts of a vector of index genes with at least the mininum combined
# score in STRING and a maximum distance (level) from the index_genes
get_all_contacts_by_level <- function(index_genes, connection_df, min_comb_score, max_level){
  # set starting level to 0
  level <- 0
  
  # create a results vector
  res <- c(index_genes)
  
  # create a vector for all genes that should be checked in the current level
  this_level_genes <- index_genes
  
  # create a vector for all genes that should be checked in the next level
  next_level_genes <- c()
  
  while (level <= max_level){
    # check all genes in the current level
    for (gene in this_level_genes){
      
      # get connection list of checked gene
      con_list <- get_direct_contacts(index_gene = gene,
                                      connection_df = connection_df,
                                      min_comb_score = min_comb_score)
      
      # fill the results vector with the current gene
      res <- union(res, gene)
      
      # add direct contacts of the current gene to next level genes (if not in results yet)
      next_level_genes <- setdiff(union(next_level_genes, con_list$direct_contacts), res)
      
      # overwrite connection df with smaller sub connection df
      connection_df <- con_list$sub_connection_df
    }
    
    # increase level
    level <- level + 1
    
    # turn next level genes to current level genes
    this_level_genes <- next_level_genes
    
    # empty next level genes for next round
    next_level_genes <- c()
  }
  return(res)
}


# function to create an edglist from all direct and indirect contacts of an index gene
# (and the index gene itself) above a minimum combined STRING score
create_edgelist <- function(connection_df, all_contacts, min_comb_score, symbol_annotation_df){
  # filter connection df for genes to plot
  interactions <- connection_df %>% 
    filter(from %in% all_contacts, to %in% all_contacts, combined_score >= min_comb_score)
  
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


# function to plot the interaction network given by an edgelist. The index gene symbols have a square marker
plot_interaction_network <- function(edgelist, disease_group_df, index_gene_symbols){
  # network() throws error if edgelist is too short => in that case, append first row of edgelist to edgelist
  if (nrow(edgelist) < 3){
    edgelist <- rbind(edgelist, edgelist[c(1), c(1:2)], edgelist[c(1), c(1:2)])
  }
  
  # create a network from the edgelist
  netw <- network(edgelist, directed = F)
  
  # get unique vertex names
  unique_vertices <- unique(c(edgelist[, 1], edgelist[, 2]))

  # count network edges
  edgecount <- network.edgecount(netw)
  
  # set edgecolor
  edgecolor <- "grey"
  
  # if no edges present, add a dummy edge (needed for adding traces to the plot later, else ERROR), set edgecolor to "white" -> not visible
  if (edgecount == 0){
    dummy_edge <- matrix(c(unique_vertices[1], "DUMMY_VERTEX"),  ncol = 2, byrow = TRUE)
    netw <- network(rbind(edgelist, dummy_edge), directed = F)
    edgecolor <- "white"
  }
  
  # add index gene info
  netw %v% "is_index_gene" = ifelse(network.vertex.names(netw) %in% index_gene_symbols, "yes", "no")
  
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
  vertex_color_map <- c(vertex_color_map, setNames("#FFFFFF", "DUMMY_VERTEX")) # color of dummy vertex is "white" -> not visible
  
  # assign colors to vertices based on their names
  vertex_colors <- sapply(network.vertex.names(netw), function(v) vertex_color_map[[v]])
  
  # create a plotly object
  p <- ggnet2(netw, 
              mode = "fruchtermanreingold",
              label = TRUE, 
              color = vertex_colors, 
              label.size = 2.5, 
              label.color = ifelse(network.vertex.names(netw) %in% unique_vertices, "black", "white"),
              edge.color = edgecolor,
              shape = "is_index_gene",
              shape.palette = c("yes" = 15, "no" = 19))
  
  p_plotly <- ggplotly(p)
  
  # remove legend for 'is_index_gene'
  p_plotly <- p_plotly %>% layout(showlegend = FALSE)
  
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
      # legendgroup = "custom_cl",
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
                              all_contacts = all_contacts_to_index, 
                              min_comb_score = min_comb_score, 
                              symbol_annotation_df = distinct(disease_group_df[c("STRING_id", "symbol")]))
  
  # plot interaction network of index gene
  interaction_plot <- plot_interaction_network(edgelist = edgelist,
                                               disease_group_df = disease_group_df,
                                               index_gene_symbols = c(index_gene))
  
  return(interaction_plot)
}


# function to create a plotly object of a network of direct and indirect contacts a of an index gene vector with
# minimum combined score in STRING and a maximum contact distance (level) of the index genes
plot_network_by_level <- function(index_genes, string_db, min_comb_score, STRING_id_vec, disease_group_df, max_level){
  
  # get all STRING interactions between all genes in STRING_id_vec with a minimum combined score in STRING
  all_interactions <- get_all_interactions_above_score(STRING_id_vec = STRING_id_vec,
                                                       string_db = string_db,
                                                       min_comb_score = min_comb_score)
  
  # get all contacts to index genes with a minimum combined score and a maximum distance
  all_contacts_by_level <- get_all_contacts_by_level(index_genes = index_genes, 
                                                     connection_df = all_interactions, 
                                                     min_comb_score = min_comb_score, 
                                                     max_level = max_level)
  
  # create an edgelist of all contacts of the index gene (direct/indirect contacts)
  edgelist <- create_edgelist(connection_df = all_interactions, 
                              all_contacts = all_contacts_by_level, 
                              min_comb_score = min_comb_score, 
                              symbol_annotation_df = distinct(disease_group_df[c("STRING_id", "symbol")]))
  
  # add edges from the index genes to themselves (as they might not have any other connections) to the edgelist
  index_gene_symbols <- disease_group_df[c("STRING_id", "symbol")] %>% 
    filter(STRING_id %in% index_genes) %>% 
    .$symbol %>% 
    unique()
  
  index_edges <- matrix(rep(index_gene_symbols, each = 2), ncol = 2, byrow = TRUE)
  edgelist <- rbind(edgelist, index_edges)
  
  # plot interaction network of index gene
  interaction_plot <- plot_interaction_network(edgelist = edgelist,
                                               disease_group_df = disease_group_df,
                                               index_gene_symbols = index_gene_symbols)
  
  return(interaction_plot)
}


# function to return x-/y-coordinates randomly drawn from a disk with given radius 'max_rad'
get_random_coordinates_on_disk <- function(max_rad){
  r <- sqrt(runif(1, min=0, max=max_rad^2) )
  alpha <- runif(1, min=0, max=2*pi)
  x_c <- r*cos(alpha)
  y_c <- r*sin(alpha)
  return(list(x=x_c, y=y_c))
}


# function to add a new point on a disk with radius 'max_rad' that has a minimum euclidean distance form the other points
add_new_point_with_min_dist <- function(coordinates_df, min_euc_dist, max_rad){
  any_true <- TRUE
  
  while(any_true){
    new_point_coord <- get_random_coordinates_on_disk(max_rad=max_rad)
    coordinates_df_m <- coordinates_df %>% 
      mutate(euc_dist = euclidean_distance(x, y, new_point_coord$x, new_point_coord$y),
             small_euc_dist = euc_dist < min_euc_dist)
    any_true <- any(coordinates_df_m$small_euc_dist)
  }
  
  coordinates_df <- rbind(coordinates_df, data.frame(x=new_point_coord$x, y=new_point_coord$y))
  return(coordinates_df)
}


# function to calculate Euclidean distance
euclidean_distance <- function(x1, y1, x2, y2) {
  sqrt((x2 - x1)^2 + (y2 - y1)^2)
}


# function to equally distribute n points on a cirlce with given radius
get_circle_coordinates <- function(n, radius=1) {
  angles <- seq(0, 2*pi, length.out = n+1)[-1]  # equally spaced angles
  coordinates <- data.frame(
    x = cos(angles) * radius,
    y = sin(angles) * radius
  )
  return(coordinates)
}


# function to get n random coordinates within a disk of radius 'max_red' with a minimum euclidean distance from each other
get_n_random_coordinates_on_disk <- function(no_cluster_genes, min_euc_dist, max_rad){
  
  # initialize the progress bar
  pb <- progress::progress_bar$new(
    format = "  progress [:bar] :percent eta: :eta",
    total = no_cluster_genes
  )
  
  # print out options in case of slow progess
  cat("\033[1;31mIf progress is too slow: Decrease 'min_euc_dist' or increase 'max_rad'!\033[0m\n")
  
  cluster_coordinates <- data.frame(x=numeric(), y=numeric())
  for (i in seq(no_cluster_genes)){
    cluster_coordinates <- add_new_point_with_min_dist(coordinates_df = cluster_coordinates, 
                                                       min_euc_dist = min_euc_dist, 
                                                       max_rad = max_rad)
    
    # increment the progress bar
    pb$tick()
  }
  
  # close the progress bar
  pb$terminate()
  
  return(cluster_coordinates)
}


# function to get coordinates of n points evenly distributed on a circle with given radius
get_cluster_center_coordinates <- function(no_clusters, cluster_center_circle_radius){
  center_coordinates <- data.frame(x=numeric(), y=numeric())
  
  if (no_clusters > 4){
    # place one cluster center in the center (0,0)
    center_coordinates <- data.frame(x=0, y=0)
    
    # equally distribute the remaining clusters in a circle with given radius
    circle_coordinates <- get_circle_coordinates(n = no_clusters - 1, radius = cluster_center_circle_radius)
    
  } else {
    # equally distribute clusters in a circle with given radius
    circle_coordinates <- get_circle_coordinates(n = no_clusters, radius = cluster_center_circle_radius)
  }
  center_coordinates <- rbind(center_coordinates, circle_coordinates)
  
  # randomly assign each cluster coordinates
  cluster_center_coordinates <- cbind(data.frame(cluster=sample(seq(no_clusters), no_clusters, replace=FALSE)), center_coordinates)
  
  return(cluster_center_coordinates)
}


# function to assign genes of a cluster to random coordinates
assign_cluster_coordinates <- function(cluster_center, cluster_genes, min_euc_dist, max_rad){
  # get coordinates of cluster genes evenly distributed on a circle with radius 'max_rad'
  cluster_coordinates <- get_n_random_coordinates_on_disk(no_cluster_genes = length(cluster_genes), 
                                                          min_euc_dist = min_euc_dist, 
                                                          max_rad = max_rad)
  
  # bind with gene symbols
  cluster_coordinates <- cbind(data.frame(gene_symbol=cluster_genes), cluster_coordinates)
  
  # add center coordinates
  cluster_coordinates <- cluster_coordinates %>% 
    mutate(x = x + cluster_center$x,
           y = y + cluster_center$y)
  
  return(cluster_coordinates)
  
}


# function to assign all genes in disease_group_df coordinates for a plot of all gene clusters
get_all_coordinates <- function(cluster_center_circle_radius,
                                min_euc_dist,
                                max_rad,
                                disease_group_df){
  
  # get main clusters
  clusters <- unique(disease_group_df$main_cluster)
  
  # get cluster center coordinates
  cluster_center_coordinates <- get_cluster_center_coordinates(no_clusters = length(clusters), 
                                                               cluster_center_circle_radius = cluster_center_circle_radius)
  
  # create empty df
  all_coordinates <- data.frame(gene_symbol = character(), x = numeric(), y = numeric())
  
  # get number of genes of largest cluster
  max_cluster_genes <- disease_group_df %>% group_by(main_cluster) %>% summarise(count=n()) %>% .$count %>% max()
  
  for (cluster in clusters){
    # get genes of this cluster
    cluster_genes <- disease_group_df %>% 
      dplyr::filter(main_cluster == cluster) %>% 
      .$symbol
    
    # get the radius of this cluster (lager cluster => larger radius, radius of largest cluster = 1)
    cluster_rad <- sqrt(length(cluster_genes)/max_cluster_genes)
    
    # get coordinates of this cluster's center
    cluster_center <- list(x = cluster_center_coordinates[cluster_center_coordinates$cluster == cluster, "x"],
                           y = cluster_center_coordinates[cluster_center_coordinates$cluster == cluster, "y"])
    
    # get coordinates for all genes of this cluster
    cluster_coordinates <- assign_cluster_coordinates(cluster_center = cluster_center, 
                                                      cluster_genes = cluster_genes,
                                                      min_euc_dist = min_euc_dist,
                                                      max_rad = cluster_rad)
    
    # add coordinates of this cluster to the dataframe with all coordinates
    all_coordinates <- rbind(all_coordinates, cluster_coordinates)
  }
  
  return(all_coordinates)
}


# function to plot all genes in clusters
plot_all_genes_in_clusters <- function(disease_group_df, cluster_center_circle_radius, min_euc_dist, max_rad){
  
  # add column with main cluster of each gene
  disease_group_df <- disease_group_df %>%
    rowwise %>%
    mutate(main_cluster = as.numeric(strsplit(cluster_index, "-")[[1]][1]))
  
  # define a custom (colorblind readable) color palette for kidney_disease_group_short
  custom_colors <- c("tubulopathy" = "#88CCEE", 
                     "glomerulopathy" = "#CC6677", 
                     "cancer" = "#DDCC77", 
                     "cakut" = "#44AA99", 
                     "cyst_cilio" = "#AA4499", 
                     "complement" = "#999933",
                     "nephrocalcinosis" = "#888888")
  
  disease_group_df$color <- custom_colors[disease_group_df$kidney_disease_group_short]
  
  all_coordinates <- get_all_coordinates(cluster_center_circle_radius = cluster_center_circle_radius,
                                         min_euc_dist = min_euc_dist,
                                         max_rad = max_rad,
                                         disease_group_df = disease_group_df)
  
  df <- all_coordinates %>% 
    left_join(disease_group_df[, c('symbol', 'color')], by=c("gene_symbol" = "symbol"), relationship = "many-to-many")
  
  # Create a scatter plot with circles and names inside
  p_plotly <- plot_ly() %>%
    add_trace(
      x = df$x,
      y = df$y,
      type = 'scatter',
      mode = 'markers',
      marker = list(symbol = 'circle', size = 30, color = df$color),
      text = df$gene_symbol,
      showlegend = FALSE,
      hoverinfo = 'text'
      
    ) %>%
    add_trace(
      x = df$x,
      y = df$y,
      type = 'scatter',
      mode = 'text',
      text = df$gene_symbol,
      textposition = 'middle center',
      textfont = list(size = 8),
      showlegend = FALSE,
      hoverinfo ='none'
      
    ) %>%
    layout(title = 'Gene clusters',
           xaxis = list(title = ''),
           yaxis = list(title = ''))
  
  
  # add legend to the plot
  p_plotly <- p_plotly %>%
    add_trace(
      x = rep(-2, length(custom_colors)),
      y = seq(-2, length.out = length(custom_colors), by = -0.15),
      type = "scatter",
      mode = "markers+text",
      marker = list(color = unname(custom_colors)), 
      text = names(custom_colors),
      marker = list(size = 12),
      showlegend = TRUE, 
      legendgroup = "vertex_colors",
      textposition = "middle right",
      textfont = list(size = 15, color = "black"),
      hoverinfo ='none'
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
  
  return(p_plotly)
}

