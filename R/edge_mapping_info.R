#' Prepare edges for mapping
#' @description Modify the mapping information for desired look when graphed in
#'  Cytoscape
#' @param expanded_edges The data frame object generated via the function 
#' expand_KEGG_edges() OR has been modified by the function add_edge_data()
#' @param data_added A logical indicator; must be set to TRUE if user data has
#'  been added (i.e. edges modified by function add_edge_data())
#' @param significance_markup A logical indicator; if set to TRUE will color 
#' edges based on direction and significance of correlation (as determined by 
#' user-data-analysis)
#' @param tidy_edge A logical indicator; must be set to FALSE for expanded edges
#' @return A data.frame object for edges that will be passed on to the function 
#' get_graph_object
#' @export
#' @examples
#' p53_KGML <- get_KGML("hsa04115")
#' p53_KEGG_mappings <- expand_KEGG_mappings(p53_KGML)
#' 
#' #Default; no data added to edges:
#' 
#' p53_edges <- expand_KEGG_edges(p53_KGML, p53_KEGG_mappings)
#' p53_edge_mapping_info <- edge_mapping_info(p53_edges) 
#' 
#' #If data is added to edges as additional attribute[s]:
#' 
#' p53_HA1E_data <- overlap_info(p53_KGML, p53_KEGG_mappings, 
#'                                "HA1E", data_type = "100_bing")
#' 
#' p53_edges_HA1E_data_MAPPED <- add_edge_data(p53_edges, p53_KEGG_mappings, 
#'                                             p53_HA1E_data, 
#'                                             data_column_no = c(3, 10,12),
#'                                             only_mapped = TRUE)
#'                                                
#' p53_edge_mapping_HA1E <- edge_mapping_info(p53_edges_HA1E_data_MAPPED, 
#'                                                       data_added = TRUE)
#'                                                       

edge_mapping_info <-
function(expanded_edges, data_added = FALSE, significance_markup = FALSE, 
         tidy_edge = TRUE){
    #expanded_edges <- expanded_edges[,-c(1)]
    for(i in 1:nrow(expanded_edges)){
        if (expanded_edges$subtype1[i] == "activation"){
            expanded_edges$color[i] = "#b20000"
        }
        else if (expanded_edges$subtype1[i] =="expression"){
            expanded_edges$color[i] = "#b20000"
        }
        else if (expanded_edges$subtype1[i] =="inhibition"){
            expanded_edges$color[i] = "#0015b2"
        } 
        else if (expanded_edges$subtype1[i] =="irreversible"){
          expanded_edges$color[i] = "#b20000"
        } 
        else if (expanded_edges$subtype1[i] =="repression"){
             expanded_edges$color[i] = "#0015b2"
        }
        else if (expanded_edges$subtype1[i] == "compound" | 
                    expanded_edges$subtype1[i] == "indirect_effect" | 
                    expanded_edges$subtype1[i] == "binding_association"){
             expanded_edges$color[i] = "#000000"
        }
        else if (expanded_edges$subtype1[i] == "de_novo") {
          expanded_edges$color[i] = "#808080"
        }
        else if (expanded_edges$subtype1[i] == "Not defined in KEGG" | 
                    expanded_edges$type[i] == "dummy") {
             expanded_edges$color[i] = "#808080"
        }
        else {
             expanded_edges$color[i] = "#1c9900"
        }
        
        if (expanded_edges$color[i] == "#1c9900"){
             expanded_edges$edge_label[i] = expanded_edges$value[i]
        }
        else {
             expanded_edges$edge_label[i] = NA
        }
    
        if (is.na(expanded_edges$subtype2[i]) == FALSE &  
            expanded_edges$subtype2[i] != "indirect"){
                if (expanded_edges$subtype1[i] == "inhibition"){
                    expanded_edges$color[i] = "#690099"
                }
                else if (expanded_edges$subtype1[i] == "activation"){
                    expanded_edges$color[i] = "#FF6600"
                }
            expanded_edges$edge_label[i] = expanded_edges$value2[i]
        }
        if(data_added){
            if (expanded_edges$premapped[i] == 1){
                if (expanded_edges$value[i] == "--|"){
                    expanded_edges$tooltip[i] <- 
                        paste0(expanded_edges$entry1symbol[i]," --| ", 
                                expanded_edges$entry2symbol[i])
                    }
                else {
                    expanded_edges$tooltip[i] <- 
                        paste0(expanded_edges$entry1symbol[i], " --> ", 
                                expanded_edges$entry2symbol[i])
                }
            }
          else {
            expanded_edges$tooltip[i] <- paste0(expanded_edges$entry1symbol[i],
                                                " - ",
                                                expanded_edges$entry2symbol[i])
            }
        }
        else {
            if (expanded_edges$value[i] == "--|"){
                expanded_edges$tooltip[i] <- 
                    paste0(expanded_edges$entry1symbol[i], " --| ",
                                                expanded_edges$entry2symbol[i])
                }
            else {
                expanded_edges$tooltip[i] <- 
                    paste0(expanded_edges$entry1symbol[i], " --> ", 
                                                expanded_edges$entry2symbol[i])
            }
        }
    }
    
    expanded_edges$name <- paste0(expanded_edges$entry1, 
                                " (", expanded_edges$type, ") ", 
                                expanded_edges$entry2)
    if (significance_markup){
        for (i in 1:nrow(expanded_edges)){
            if (expanded_edges$has_data[i] == 1){
                if (expanded_edges$summary_score[i] > 0 & 
                    expanded_edges$significant[i] == 1){
                        expanded_edges$color[i] <- "#b20000"
                        }
                else if (expanded_edges$summary_score[i] > 0 & 
                        expanded_edges$significant[i] == 0){
                        expanded_edges$color[i] <- "#FF6600"
                        }
                else if (expanded_edges$summary_score[i] <= 0 & 
                        expanded_edges$significant[i] == 1){
                        expanded_edges$color[i] <- "#0015b2"
                        }
                else {
                        expanded_edges$color[i] <- "#690099"
                }
            }
            else {
                expanded_edges$color[i] = "#808080"
            }
        }
    }
    #edge_map <- edge_map[,c(2,5,3,1,4:ncol(edge_map))]
    expanded_edges<- cbind(expanded_edges[,c("entry1", "entry2","edgeID", 
                                        "entry1accession", "entry2accession")],
                          expanded_edges[,c(6:ncol(expanded_edges))])
    expanded_edges <- expanded_edges[order(expanded_edges$edgeID),]

    if (tidy_edge == TRUE) {
      edge_IDs <- seq(min(expanded_edges$edgeID), max(expanded_edges$edgeID))
      for (i in edge_IDs[edge_IDs %in% expanded_edges$edgeID]){
        if(data_added == TRUE){
          expanded_edges <- tidy_edge(edges = expanded_edges,
                                      edge_id = edge_IDs[i], 
                                      data_added = TRUE,
                                      by_significance = TRUE)
        }
        if(data_added == FALSE){
          expanded_edges <- tidy_edge(edges = expanded_edges,
                                      edge_id = edge_IDs[i], 
                                      data_added = FALSE)
        }
      }
    }
    return(expanded_edges)
}
