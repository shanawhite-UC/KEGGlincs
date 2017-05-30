#' Tidy up pathway by combining edges inside of edge_mapping_info
#' @description Combine edges that share nodes and have other commonalities
#' @export
#' @param edges The edge dataframe 
#' @param edge_id The numeric value for the edge_id
#' @param by_significance A logical indicator; option if data is added
#' @param by_number A logical indicator; gives rough estimate of edge amount
#' @return A data frame that has had the given edge condensed for viewing
#' @examples \dontrun{
#' if (tidy_edge == TRUE) {
#'    edge_IDs <- seq(min(expanded_edges$edgeID), max(expanded_edges$edgeID))
#'    for (i in edge_IDs){
#'      if(data_added == TRUE){
#'        expanded_edges <- tidy_edge(edges = expanded_edges,
#'                                    edge_id = edge_IDs[i], 
#'                                    data_added = TRUE,
#'                                    by_significance = TRUE)
#'      }
#'      if(data_added == FALSE){
#'        expanded_edges <- tidy_edge(edges = expanded_edges,
#'                                    edge_id = edge_IDs[i], 
#'                                    data_added = FALSE)
#'      }
#'    }
#'
#'}
#'}

tidy_edge <- function(edges, edge_id, data_added = TRUE, 
                      by_significance = FALSE, by_number = TRUE){
  edge <- edges[edges$edgeID == edge_id,]
  edges_1 <- edges[edges$edgeID != edge_id,]
  if (!data_added){
    entry1accessions <- paste(unique(edge$entry1accession), collapse =",")
    entry2accessions <- paste(unique(edge$entry2accession), collapse =",")
    entry1symbols <- paste(unique(edge$entry1symbol), collapse = ",")
    entry2symbols <- paste(unique(edge$entry2symbol), collapse = ",")
    tooltip <- paste(entry1symbols, edge$value[1], entry2symbols)
    edge$entry1symbol[1] <- entry1symbols
    edge$entry2symbol[1] <- entry2symbols
    edge$entry1accession[1] <- entry1accessions
    edge$entry2accession[1] <- entry2accessions
    edge$tooltip[1] <- tooltip
    reduced_edge <- edge[1,]
  }
  if (data_added){
    if (sum(edge$has_data) == 0){
      entry1accessions <- paste(unique(edge$entry1accession), collapse =",")
      entry2accessions <- paste(unique(edge$entry2accession), collapse =",")
      entry1symbols <- paste(unique(edge$entry1symbol), collapse = ",")
      entry2symbols <- paste(unique(edge$entry2symbol), collapse = ",")
      tooltip <- paste(entry1symbols, edge$value[1], entry2symbols)
      edge$entry1symbol[1] <- entry1symbols
      edge$entry2symbol[1] <- entry2symbols
      edge$entry1accession[1] <- entry1accessions
      edge$entry2accession[1] <- entry2accessions
      edge$tooltip[1] <- tooltip
      reduced_edge <- edge[1,]
    }
    
    if (sum(edge$has_data) != 0 ){
      reduced_edge <- edge[edge$has_data != 0,]
      if(by_significance == TRUE){
        for (i in 1:length(unique(reduced_edge$color))){
          rp <- reduced_edge[reduced_edge$color == unique(reduced_edge$color)[i],]
          reduced_edge <- reduced_edge[reduced_edge$color != unique(reduced_edge$color)[i],]
          entry1accessions <- paste(unique(rp$entry1accession), collapse =",")
          entry2accessions <- paste(unique(rp$entry2accession), collapse =",")
          entry1symbols <- paste(unique(rp$entry1symbol), collapse = ",")
          entry2symbols <- paste(unique(rp$entry2symbol), collapse = ",")
          tooltip <- paste(unique(rp$tooltip), collapse = " , ")
          average_summary_score <- mean(rp$summary_score)
          rp$entry1symbol[1] <- entry1symbols
          rp$entry2symbol[1] <- entry2symbols
          rp$entry1accession[1] <- entry1accessions
          rp$entry2accession[1] <- entry2accessions
          rp$tooltip[1] <- tooltip
          rp$summary_score[1] <- average_summary_score
          
          if(by_number == TRUE) {
            col <- rp$color[1]
            thresh <- min(5, (nrow(rp)-1))
            value <- 1-(0.1*thresh)
            col1 <- round(value*(col2rgb(col)[1]))
            col2 <- round(value*(col2rgb(col)[2]))
            col3 <- round(value*(col2rgb(col)[3]))
            rp$color[1] <- rgb(col1, col2, col3, maxColorValue = 255)
          }    
          reduced_edge <- rbind(rp[1,], reduced_edge)
        }
      }
    }
  }
  edges <- rbind (edges_1, reduced_edge)
  return(edges)
}