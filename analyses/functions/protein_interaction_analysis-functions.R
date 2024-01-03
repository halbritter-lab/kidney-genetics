
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