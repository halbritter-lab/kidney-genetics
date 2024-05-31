# TODO: check which libraries to omit
# TODO: combine with D_ProteinInteractionAnalysis.R
# TODO: clean
# TODO: change plot_interaction_network(), so that the input is already a network object
# TODO: move functions to function file

library(plotly)
library(ggnet)
library(network)
library(sna)
library(ggplot2)
library(GGally)

string_db_full <- STRINGdb::STRINGdb$new(version = string_db_version,
                                         species = 9606,
                                         score_threshold = 100,
                                         input_directory = string_db_files_path)

# get STRING IDs of genes that have a kidney disease group
STRING_id_vec <- unique(kid_groups$STRING_id)

min_comb_score <- 400
index_gene <- "9606.ENSP00000215832"

all_interactions <- string_db_full$get_interactions(STRING_id_vec) %>% 
  distinct() %>% 
  filter(combined_score >= min_comb_score)

# function that checks if two genes are connected with a minimum combined score: works, but slow
are_connected <- function(gene1, gene2, min_comb_score){
  con <- string_db_full$get_interactions(c(gene1, gene2)) %>% 
    filter(combined_score >= min_comb_score) %>% 
    nrow()
  
  if (con > 0){return(1)}
  else{return(0)}
}

# other version - faster
are_connected <- function(gene1, gene2, min_comb_score){
  df <- all_interactions
  df$correct_connection <- ifelse((df$from == gene1 & df$to == gene2 & df$combined_score >= min_comb_score) |
                                     (df$to == gene1 & df$from == gene2 & df$combined_score >= min_comb_score),
                                   1, 0)

  if (1 %in% df$correct_connection){return(1)}
  else{return(0)}
}

# examples
are_connected("9606.ENSP00000215832", "9606.ENSP00000219476", min_comb_score = 800)
are_connected("9606.ENSP00000219476", "9606.ENSP00000215832", min_comb_score = 960)

get_direct_contacts <- function(index_gene, min_comb_score){
  rem_genes <- setdiff(STRING_id_vec, index_gene)
  connected_genes <- c()
  for (i in rem_genes){
    if (are_connected(index_gene, i, min_comb_score)){
      connected_genes <- c(connected_genes, i)
    }
  }
  return(connected_genes)
}

get_direct_contacts(index_gene = "9606.ENSP00000215832", min_comb_score = 400)

# get all direct contacts and indirect contacts of an index gene: works, but slow!!!
get_all_contacts <- function(index_gene, min_comb_score, res = c()){
  # add index gene to results
  res <- union(res, index_gene)
  
  # get directly connected genes
  direct_contacts <- get_direct_contacts(index_gene, min_comb_score)
  
  # get genes that have not been checked for directly connected genes
  newly_added <- setdiff(direct_contacts, res)
  
  # if (length(res) > 50){c(res, direct_contacts)}
  
  if (length(newly_added) == 0){
    return(c(res))
    } else {
      for (i in newly_added){
        # res <- c(res, i)
        print(length(res))
        res <- get_all_contacts(i, min_comb_score, res)
        # res <- union(res, Recall(i, min_comb_score, res = res))
      }
      return(res)
  }
}
    
# example
res1 <- get_all_contacts(index_gene = "9606.ENSP00000215832", min_comb_score = 900, res=c())

# other version without recursive call
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


all_interactions <- all_interactions %>% 
  mutate(both_genes = list(c(from, to))) %>% 
  mutate(row_index = cur_group_id())


# get direct contacts of a gene and remove corresponding rows in df
get_direct_contacts_v2 <- function(index_gene, connection_df, min_comb_score){
  gene_connections <- connection_df %>% filter(index_gene %in% both_genes, combined_score >= min_comb_score)
  direct_contacts <- setdiff(union(gene_connections$from, gene_connections$to), index_gene)
  # rows_to_remove <- connection_df$row_index
  sub_connection_df <- connection_df %>% filter(!(row_index %in% gene_connections$row_index))
  return(list(direct_contacts = direct_contacts, sub_connection_df = sub_connection_df))
}


index_gene <- "9606.ENSP00000290866"
connection_df <- all_interactions
min_comb_score <- 850

# get all direct and indirect contacts of an index gene (v2: without recursive call, probably faster)
get_all_contacts_v2 <- function(index_gene, connection_df, min_comb_score){
  # get direct contacts of index_gene
  con_list <- get_direct_contacts_v2(index_gene, connection_df, min_comb_score)
  
  # create a results vector
  res <- c(index_gene)
  
  # create a vector for genes that have to be checked for further contacts
  genes_to_check <- con_list$direct_contacts
  
  # check all genes in genes_to_check for further contacts
  while (length(genes_to_check) > 0){
    gene1 <- genes_to_check[1]
    
    # get direct contacts of gene that is checked
    con_list <- get_direct_contacts_v2(gene1, con_list$sub_connection_df, min_comb_score)
    
    # add gene that is checked to results vector
    res <- union(res, gene1)
    
    # removed checked gene from genes_to_check
    genes_to_check <- setdiff(genes_to_check, gene1)
    
    # add new genes (direct contacts of checked gene) to genes_to_check
    genes_to_check <- union(genes_to_check, setdiff(con_list$direct_contacts, res))
  }
  
  return(res)
}


# example
res <- get_all_contacts_v2(index_gene = "9606.ENSP00000290866",  # "9606.ENSP00000355627"

                           connection_df = all_interactions,
                           min_comb_score = 850)

min_comb_score <- 850
#############################

# TODO: add maxgroups_full, hgnc_annotated

plot_interaction_network <- function(all_interaction_df, STRING_id_vec, min_comb_score){
  # get interactions for STRING_id_vec
  netw <- all_interaction_df %>% filter(from %in% res, to %in% res, combined_score > min_comb_score)
  
  # get symbols
  netw <- netw %>% 
    left_join(hgnc_annotated[c('STRING_id', 'symbol')], by = c("from" = "STRING_id")) %>% rename(symbol_from = symbol) %>% 
    left_join(hgnc_annotated[c('STRING_id', 'symbol')], by = c("to" = "STRING_id")) %>% rename(symbol_to = symbol) %>% distinct() %>% 
    dplyr::select(symbol_from, symbol_to, combined_score)  
  
  # create an edgelist
  edgel <- cbind(netw$symbol_from, netw$symbol_to)
  
  # create a network object
  net2 <- network(edgel, directed = F)
  
  # get unique vertex names
  unique_vertices <- unique(c(edgel[, 1], edgel[, 2]))
  
  
  # define a custom (colorblind readable) color palette for kidney_disease_group_short
  custom_colors <- c("tubulopathy" = "#88CCEE", 
                     "glomerulopathy" = "#CC6677", 
                     "cancer" = "#DDCC77", 
                     "cakut" = "#44AA99", 
                     "cyst_cilio" = "#AA4499", 
                     "complement" = "#999933",
                     "nephrocalcinosis" = "#888888")
  
  # define vertex colors
  vertex_colors <- max_groups_full %>% 
    filter(STRING_id %in% res) %>% 
    mutate(color = custom_colors[kidney_disease_group_short])
  
  vertex_color_map <- setNames(vertex_colors$color, vertex_colors$symbol)
  
  # assign colors to vertices based on their names
  vertex_colors <- sapply(network.vertex.names(net2), function(v) vertex_color_map[[v]])
  
  # create a plotly object
  p <- ggnet2(net2, label = TRUE, color = vertex_colors, label.size = 2.5)
  
  p_plotly <- ggplotly(p)
  
  # add legend to the plot
  p_plotly <- p_plotly %>%
    add_trace(
      x = rep(-0.1, length(custom_colors)),
      y = seq(-0.2, length.out = length(custom_colors), by = -0.1),
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


# example
junk <- plot_interaction_network(all_interaction_df = all_interactions,
                         STRING_id_vec = res,
                         min_comb_score = 850
                         )

junk

# NOTES - TO BE DELETED
# junk <- string_db_full$get_interactions(res) %>% filter(combined_score >= 850)
# 
# union(junk$from, junk$to) %>% unique %>% length


# # get also reverse interactions so that the interaction matrix later will be symmetrical
# netw_copy <- netw %>% rename(from = to, to = from)
# netw <- rbind(netw, netw_copy)




# 
# # Create a trace for the legend points
# legend_trace <- list(
#   x = rep(1.1, 7),
#   y = c(1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7),
#   mode = "markers",
#   marker = list(color = c("#FF0000", "#FF8B00", "#E8FF00", "#5DFF00", "#00FF2E", "#00FFB9", "#00B9FF")),
#   text = names(custom_colors),
#   showlegend = TRUE, 
#   legendgroup = "vertex_colors",
#   name = "Vertex Colors"
# )
# 
# # Add the legend trace and update the layout
# p_plotly <- p_plotly %>%
#   add_trace(legend_trace) %>%
#   layout(
#     showlegend = TRUE
#   )
# 
# 
# p_plotly
# 

# Assuming you've already defined the plot 'p' using ggnet2 and converted it to plotly 'p_plotly'

# # Create a custom legend with shapes
# legend_data <- data.frame(
#   x = rep(1.1, 7), # Adjust x position
#   y = 1:7
# )



# Assuming 'p_plotly' contains your plotly object
# Assuming 'p_plotly' contains your plotly object
# 
# # Define marker colors and names
# marker_colors <- c("#FF0000", "#FF8B00", "#E8FF00", "#5DFF00", "#00FF2E", "#00FFB9", "#00B9FF")
# marker_names <- names(custom_colors)
# 
# # Create a trace for the legend points
# legend_trace <- list(
#   x = rep(1.1, length(marker_names)),
#   y = seq(1.1, 1.1 + length(marker_names) - 1, by = 0.1),  # Adjust y positions
#   mode = "markers+text",
#   marker = list(color = marker_colors, size = 10),
#   text = marker_names,
#   textposition = "bottom center",
#   showlegend = TRUE,
#   legendgroup = "vertex_colors",
#   name = "Vertex Colors"
# )
# 
# # Add the legend trace and update the layout
# p_plotly <- p_plotly %>%
#   add_trace(legend_trace) %>%
#   layout(
#     showlegend = TRUE
#   )
# 

# custom_colors <- c("tubulopathy" = "#FF0000", 
#                    "glomerulopathy" = "#FF8B00", 
#                    "cancer" = "#E8FF00", 
#                    "cakut" = "#5DFF00", 
#                    "cyst_cilio" = "#00FF2E", 
#                    "complement" = "#00FFB9",
#                    "nephrocalcinosis" = "#00B9FF")  


# p_plotly



# # Recursive function to get all contacts
# get_all_contacts <- function(index_person, min_comb_score, visited = c()) {
#   # Add the current index_person to visited
#   visited <- c(visited, index_person)
#   
#   # Get direct contacts of the current index_person
#   direct_contacts <- get_direct_contacts(index_person, min_comb_score)
#   
#   # Find contacts that have not been visited
#   new_contacts <- setdiff(direct_contacts, visited)
#   
#   # Base case: If there are no new contacts, return the visited list
#   if (length(new_contacts) == 0) {
#     return(visited)
#   } else {
#     # Recursively get contacts for each new contact
#     for (new_contact in new_contacts) {
#       print(length(visited))
#       visited <- get_all_contacts(new_contact, min_comb_score, visited)
#     }
#     return(visited)
#   }
# }

# res1 <- get_all_contacts(index_person = "9606.ENSP00000215832", min_comb_score = 900, visited=c())


