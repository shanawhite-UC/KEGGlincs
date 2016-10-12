#' Combines all other package functions for one-step pathway visualization
#' @export
#' @param pathwayid A KEGG pathway ID of the form "hsa12345" 
#' (only human pathways currently)
#' @param cell_line If left as NA will generate a pathway map without 
#' data-dependent attributes (such as edge width).  To use in combination 
#' with LINCS data, choose from the set of cell lines: 
#' (A375,A549,ASC,HA1E,HCC515,HEK293T,HEKTE,HEPG2,HT29,MCF7,NCIH716,NPC,PC3,
#' SHSY5Y,SKL,SW480,VCAP)
#' @param refine_by_cell_line A logical indicator
#' @param add_L1000_edge_data A logical indicator 
#' @param data_type Choose from data types: (100_full, 100_bing, 50_lm)
#' @param pert_time Choose from (6,24,48,96,120,144,168)
#' @param only_mapped A logical indicator; if set to FALSE will return 'de-novo'
#' edges that 'exist' in data but are not documented in KEGG
#' @param significance_markup A logical indicator; if set to TRUE will color
#'  edges based on direction and significance of correlation (as determined by 
#'  user-data-analysis)
#' @param layered_nodes A logical indicator; if set to TRUE will create a graph 
#' with 'stacked' nodes that the user can manipulate when multiple nodes are 
#' mapped to one location
#' @param graph_title An optional user-specified graph title
#' @param get_data A logical indicator; if set to true, will return the 
#' 'expanded' edge information for the specified pathway
#' @param convert_KEGG_IDs A logical indicator; if set to TRUE KEGG 
#' compounds will remain labeled via KEGG codes (do not need KEGGREST)
#' @return  A dynamic map in Cytoscape automatically formatted for convenient 
#' viewing and, if idicated by user, a data.frame object with detailed 
#' information for 'expanded' KEGG edges
#' @examples \dontrun{ 
#' 
#' #Default KEGG pathway with colored edges representing type of relationship:
#' KEGG_lincs("hsa04115", convert_KEGG_IDs = FALSE)
#' 
#' #KEGG pathway with edge width and color based on observed experimental data:
#' KEGG_lincs("hsa04115", "HA1E")
#' 
#' #Have edge information data.frame to be output to the global environment:
#' p53_edge_info <- KEGG_lincs("hsa04115", graph_title = "p53"
#'                              convert_KEGG_IDs = FALSE, get_data = TRUE)
#'                              }


KEGG_lincs <-
function(pathwayid, cell_line = NA,
                            refine_by_cell_line = NA,
                            add_L1000_edge_data = TRUE,  
                            significance_markup = TRUE,
                            data_type = "100_full",
                            pert_time = 96,
                            only_mapped = TRUE,
                            layered_nodes = FALSE, graph_title = "default",
                            get_data = FALSE,
                            convert_KEGG_IDs = TRUE){
    if (is.na(refine_by_cell_line)){
        if (is.na(cell_line)){
            refine_by_cell_line <- FALSE
        }
        else {
            refine_by_cell_line <- TRUE
        }
    }

    KGML <- get_KGML(pathwayid, get_if_no_edges = TRUE)
    if(!isS4(KGML)){
        return()
    }
    KEGG_mappings <- 
        expand_KEGG_mappings(KGML,convert_KEGG_IDs = convert_KEGG_IDs)
    if (refine_by_cell_line) {
        full_mappings <- KEGG_mappings
        KEGG_mappings <- refine_mappings(KEGG_mappings, cell_line)
        for (i in 1:nrow(full_mappings)){
            if(!full_mappings$entryID[i] %in% KEGG_mappings$entryID){
                full_mappings$BGcolor[i] <- "#d3d3d3"
                full_mappings$in_relationship <- 0
            }
        }
    }
    expanded_edges <- expand_KEGG_edges(KGML, KEGG_mappings)
    if(expanded_edges$type[1] == "dummy"){
        graph_title <- paste0("Pathway = ", pathwayid, ":", 
                        KGML@pathwayInfo@title, "Cell-Line: ", cell_line, 
                        "  *No Edges in Pathway")
    }

    if (is.na(cell_line)){
        edge_map <- edge_mapping_info(expanded_edges)
    if (graph_title == "default"){
        graph_title <- paste0("Pathway = ", pathwayid, ":", 
                            KGML@pathwayInfo@title)
        }
    }
    if (!is.na(cell_line) & !add_L1000_edge_data){
        edge_map <- edge_mapping_info(expanded_edges)
        if (graph_title == "default"){
            graph_title <- paste0("Pathway = ", pathwayid, ":", 
                            KGML@pathwayInfo@title, "Refined by Cell-Line:  ",
                            cell_line)
        }
    }
    if (nrow(expanded_edges[expanded_edges$type == "maplink",]) == 
        nrow(expanded_edges) & only_mapped) {
            edge_map <- edge_mapping_info(expanded_edges)
            if (graph_title == "default"){
                graph_title <- paste0("Pathway = ", pathwayid, ":", 
                            KGML@pathwayInfo@title, "Cell-Line: ", 
                            cell_line,  "  *No Edges in Data")
            }
        warning("All documented edges are of type 'maplink'; 
            Overlap data cannot be mapped to selected pathway")
    }
    if(!is.na(cell_line) & add_L1000_edge_data) {
        user_data <- overlap_info(KGML, KEGG_mappings, cell_line, 
                                data_type = data_type, pert_time = pert_time)
        if (!is.na(user_data)[1,1]) {
            edges_plus_data <- add_edge_data(expanded_edges, KEGG_mappings, 
                                            user_data, c(10,12), 
                                            only_mapped = only_mapped)
            edge_map <- edge_mapping_info(edges_plus_data, data_added = TRUE, 
                                    significance_markup = significance_markup)
            if (graph_title == "default"){
                graph_title <- paste0("Pathway = ", pathwayid, ":", 
                                KGML@pathwayInfo@title, ",  Cell-Line: ", 
                                cell_line, ",  Data type: ", data_type)
            }
        }
        else {
            edge_map <- edge_mapping_info(expanded_edges)
            if (graph_title == "default"){
                graph_title <- paste0("Pathway = ", pathwayid, ":", 
                                KGML@pathwayInfo@title, ",  Cell-Line: ", 
                                cell_line, "   *No Edges in Data")
            }
        }
    }
    if(refine_by_cell_line){
        node_map <- node_mapping_info(full_mappings)
    }
    else{
        node_map <- node_mapping_info(KEGG_mappings)
    }
    graph_object <- get_graph_object(node_map, edge_map, 
                                    layered_nodes = layered_nodes)
    cyto_vis(graph_object, title = graph_title)
    if(get_data){
        return(edge_map)
    }
}
