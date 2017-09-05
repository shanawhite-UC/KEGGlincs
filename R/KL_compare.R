#' Combines all other package functions for one-step cell line comparison
#' @export
#' @param pathwayid A KEGG pathway ID of the form "hsa12345" 
#' (only human pathways currently)
#' @param cell_line1 Choose from the set of cell lines: 
#' (A375,A549,ASC,HA1E,HCC515,HEK293T,HEKTE,HEPG2,HT29,MCF7,NCIH716,NPC,PC3,
#' SHSY5Y,SKL,SW480,VCAP)
#' @param cell_line2 A cell line such that cell_line1 != cell_line2
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
#' # Compare p53 pathway between cell lines A375 and A549:
#' KL_compare("hsa04115", "A375", "A549")
#'}

KL_compare <-
    function(pathwayid, cell_line1 = NA, cell_line2 = NA,
             refine_by_cell_line = TRUE,
             add_L1000_edge_data = TRUE,  
             significance_markup = TRUE,
             data_type = "100_full",
             pert_time = 96,
             only_mapped = TRUE,
             get_data = FALSE,
             convert_KEGG_IDs = TRUE,
             graph_title = "default",
             tidy_edges = TRUE,
             layered_nodes = FALSE){
        cell_lines <- c("A375","A549","ASC","HA1E","HCC515","HEK293T","HEKTE",
                        "HEPG2","HT29","MCF7","NCIH716","NPC","PC3","SHSY5Y",
                        "SKL","SW480","VCAP")
        if (is.na(cell_line1)){
            warning(paste0(
                'Please choose one of the following for cell_line1: ', 
                cell_lines))
            return()
        }
        if (is.na(cell_line2)){
            warning(paste0('Please choose one the following for cell_line2: ',
                           list(cell_lines[cell_lines != cell_line1])))
            return()
        }
        
        if (cell_line2 == cell_line1){
            warning(paste0('cell_line1 = cell_line2; ', 
                           'please choose from the following for cell_line2: ',
                           list(cell_lines[cell_lines != cell_line1])))
            return()
        }
        
        KGML <- get_KGML(pathwayid, get_if_no_edges = TRUE)
        if(!isS4(KGML)){
            return()
        }
        KEGG_mappings <- 
            expand_KEGG_mappings(KGML,convert_KEGG_IDs = convert_KEGG_IDs)
        if (refine_by_cell_line) {
            full_mappings1 <- KEGG_mappings
            KEGG_mappings1 <- refine_mappings(KEGG_mappings, cell_line1)
            for (i in 1:nrow(full_mappings1)){
                if(!full_mappings1$entryID[i] %in% KEGG_mappings$entryID){
                    full_mappings1$BGcolor[i] <- "#d3d3d3"
                    full_mappings1$in_relationship <- 0
                }
            }
            full_mappings2 <- KEGG_mappings
            KEGG_mappings2 <- refine_mappings(KEGG_mappings, cell_line2)
            for (i in 1:nrow(full_mappings2)){
                if(!full_mappings2$entryID[i] %in% KEGG_mappings$entryID){
                    full_mappings2$BGcolor[i] <- "#d3d3d3"
                    full_mappings2$in_relationship <- 0
                }
            }
        }
        expanded_edges <- expand_KEGG_edges(KGML, KEGG_mappings)
        if(expanded_edges$type[1] == "dummy"){
            graph_title <- paste0("Pathway = ", pathwayid, ":", 
                                  KGML@pathwayInfo@title,
                                  "Cell-Line: ", cell_line, 
                                  "  *No Edges in Pathway")
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
        cell_data1 <- overlap_info(KGML, KEGG_mappings, cell_line1, 
                                   data_type = data_type, 
                                   pert_time = pert_time)
        
        cell_data2 <- overlap_info(KGML, KEGG_mappings, cell_line2, 
                                   data_type = data_type, 
                                   pert_time = pert_time)
        
        for (i in 1:nrow(cell_data1)){
            UP <- cell_data1$UP[i] + 0.5
            DOWN <- cell_data1$DOWN[i] + 0.5
            UK1_DK2 <- cell_data1$UK1_DK2[i] + 0.5
            DK1_UK2 <- cell_data1$DK1_UK2[i] + 0.5
            cell_data1$OR_1[i] <- (UP*DOWN)/(UK1_DK2*DK1_UK2)
            cell_data1$SE_1[i] <-sqrt(1/UP + 1/DOWN + 1/UK1_DK2 + 1/DK1_UK2)
            cell_data1$log_OR_1[i] <- log(cell_data1$OR_1[i])
        }
        
        for (i in 1:nrow(cell_data2)){
            UP <- cell_data2$UP[i] + 0.5
            DOWN <- cell_data2$DOWN[i] + 0.5
            UK1_DK2 <- cell_data2$UK1_DK2[i] + 0.5
            DK1_UK2 <- cell_data2$DK1_UK2[i] + 0.5
            cell_data2$OR_2[i] <- (UP*DOWN)/(UK1_DK2*DK1_UK2)
            cell_data2$SE_2[i] <-sqrt(1/UP + 1/DOWN + 1/UK1_DK2 + 1/DK1_UK2)
            cell_data2$log_OR_2[i] <- log(cell_data2$OR_2[i])
        }
        
        if (!is.na(cell_data1)[1,1] & !is.na(cell_data2)[1,1]) {
            edges_plus_data1 <- add_edge_data(expanded_edges, KEGG_mappings, 
                                              cell_data1, c(15,16), 
                                              only_mapped = only_mapped)
            edges_plus_data2 <- add_edge_data(expanded_edges, KEGG_mappings, 
                                              cell_data2, c(15,16), 
                                              only_mapped = only_mapped)
            edges_plus_data1$unique_ID <- paste0(edges_plus_data1$entry1symbol, ":", edges_plus_data1$entry2symbol, ":", edges_plus_data1$edgeID)
            edges_plus_data2$unique_ID <- paste0(edges_plus_data2$entry1symbol, ":", edges_plus_data2$entry2symbol, ":", edges_plus_data2$edgeID)
            
            edges_compare <- merge(edges_plus_data1, edges_plus_data2[,c(18,19,21)], by = "unique_ID")
            
            for(i in 1:nrow(edges_compare)){
                if (!is.na(edges_compare$log_OR_1[i]) & !is.na(edges_compare$log_OR_2[i])){
                    edges_compare$test[i] <- (edges_compare$log_OR_1[i] - edges_compare$log_OR_2[i])/
                        sqrt(edges_compare$SE_1[i]^2 + edges_compare$SE_2[i]^2)
                    edges_compare$summary_score[i] <- exp(abs(edges_compare$test[i]))
                    
                    
                    if (edges_compare$test[i] <= qnorm(0.1) | edges_compare$test[i] >= qnorm(0.9)){
                        edges_compare$significant[i] <- 1
                    }
                    else {
                        edges_compare$significant[i] <- 0
                    }
                }
                else {
                    edges_compare$test[i] <- NA
                    edges_compare$summary_score[i] <- 0
                    edges_compare$significant[i] <- NA
                }
            }
            
            edge_map <- edge_mapping_info(edges_compare, data_added = TRUE)
            for (i in 1:nrow(edge_map)){
                if (is.na(edge_map$test[i])){
                    edge_map$has_data[i] <- 0
                    edge_map$color[i] <- "#808080"
                }
            }
            if ("premapped" %in% names (edge_map)) {
                premapped <- edge_map$premapped
                drop <- "premapped"
                edge_map <- edge_map[, ! names(edge_map) %in% drop]
                edge_map$premapped <- premapped
            }
            for (i in 1:nrow(edge_map)){
                if (!is.na(edge_map$log_OR_1[i]) & 
                    !is.na(edge_map$log_OR_2[i])){
                    if (edge_map$test[i] < 0){
                        if (edge_map$significant[i] == 0){
                            edge_map$color[i] <- "#82E0AA"
                        }
                        if (edge_map$significant[i] == 1){
                            edge_map$color[i] <- "#1E8449"
                        }
                    }
                    if (edge_map$test[i] > 0){
                        if (edge_map$significant[i] == 0){
                            edge_map$color[i] <- "#E67E22"
                        }
                        if (edge_map$significant[i] == 1){
                            edge_map$color[i] <- "#BA4A00"
                        }
                    }
                }
            }
            if (tidy_edges == TRUE) {
                edge_IDs <- seq(min(edge_map$edgeID), max(edge_map$edgeID))
                for (i in edge_IDs){
                    edge_map <- tidy_edge(edges = edge_map, 
                                          edge_id = edge_IDs[i],
                                          by_significance = TRUE)
                }
            }
            
            if (graph_title == "default"){
                graph_title <- paste0("Pathway = ", pathwayid, ":",
                                      KGML@pathwayInfo@title, " - ",
                                      cell_line1 ,"vs", cell_line2, ", 
                                      Data type: ", data_type)
            }
            }
        
        node_map <- node_mapping_info(KEGG_mappings)
        
        graph_object <- get_graph_object(node_map, edge_map, 
                                         layered_nodes = layered_nodes)
        
        edge_width_attribute = "summary_score"
        
        if (edge_width_attribute %in% names(igraph::edge_attr(graph_object))){
            min.summary_score <- min(abs(igraph::E(graph_object)$summary_score),
                                     na.rm = TRUE)
            max.summary_score <- max(abs(igraph::E(graph_object)$summary_score), 
                                     na.rm = TRUE)
            map_edge_width <- TRUE
        }
        
        cyto_vis(graph_object, title = graph_title, 
                 edge_width_attribute = "summary_score")
        if(get_data){
            return(edge_map)
        }
    }
