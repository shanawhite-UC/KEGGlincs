#' Add in edges to map documented in other pathways
#' @description For a specific pathway entity(gene), search KEGG databases
#' to see if it has any other documented relationships in KEGG.
#' expand_KEGG_edges 
#' @export
#' @param entry_accession The Accession # of the pathway entity to 'keggerize'
#' @param KGML The KGML file of the current pathway
#' @param KEGG_mappings KEGG mappings for the current pathway
#' @param edges The expanded edges for the current pathway
#' @return A modified expanded edges data frame with additional rows for new 
#' entries
#' @examples \dontrun{
#' KGML <- get_KGML("hsa04150")
#' KEGG_mappings <- expand_KEGG_mappings(KGML)
#' edges <- expand_KEGG_edges(KGML, KEGG_mappings)
#' entry_accession <- "2475"
#' mtor_plus_mtor <- keggerize_edges(entry_accession = entry_accession, 
#'                                   KGML = KGML,KEGG_mappings = KEGG_mappings,
#'                                   edges = edges) }
#' 

keggerize_edges <- function(entry_accession, KGML, KEGG_mappings, edges){
    if(!"primary_pathway" %in% names(edges)){
        edges$primary_pathwayID <- strsplit(KGML@pathwayInfo@name, ":")[[1]][2]
        edges$primary_pathway <- KGML@pathwayInfo@title
    }
    entry_code <- paste0("hsa:",entry_accession)
    entry_info <- keggGet(entry_code)
    entry_pathways <- entry_info[[1]]["PATHWAY"]
    entry_pathways <- names(entry_pathways$PATHWAY[1:length(unlist(entry_pathways))])
    current_pathway <- strsplit(KGML@pathwayInfo@name, ":")[[1]][2]
    entry_pathways <- entry_pathways[!entry_pathways %in% current_pathway]
    
    for(p in 1:length(entry_pathways)){
        pathway_KGML <- get_KGML(entry_pathways[p])
        pathway_KEGG_mappings <- expand_KEGG_mappings(pathway_KGML)
        pathway_edges <- expand_KEGG_edges(pathway_KGML, pathway_KEGG_mappings)
        entry_edges <- pathway_edges[pathway_edges$entry1accession == entry_accession |
                                         pathway_edges$entry2accession == entry_accession, ]
        entry_sources <- entry_edges[entry_edges$entry1accession != entry_accession &
                                         entry_edges$entry1accession %in% KEGG_mappings$entryACCESSION,]
        entry_sources <- entry_sources[!duplicated(entry_sources[,6:16]),]
        if(nrow(entry_sources) > 0){
            for(i in 1:nrow(entry_sources)){
                entry1_locations <- KEGG_mappings$entryID[KEGG_mappings$entrySYMBOL == entry_sources$entry1symbol[i]]
                entry2_locations <- KEGG_mappings$entryID[KEGG_mappings$entrySYMBOL == entry_sources$entry2symbol[i]]
                if(length(entry1_locations) > 1 | length(entry2_locations > 1)){
                    
                    distance_matrix <- matrix(nrow = length(entry2_locations),
                                              ncol = length(entry1_locations))
                    colnames(distance_matrix) <- entry1_locations
                    rownames(distance_matrix) <- entry2_locations
                    for(j in 1:nrow(distance_matrix)){
                        for (k in 1:ncol(distance_matrix)){
                            value <-
                                sqrt((KEGG_mappings$Xcoord[KEGG_mappings$entryID ==
                                                               colnames(distance_matrix)[k]][1]
                                      - KEGG_mappings$Xcoord[KEGG_mappings$entryID ==
                                                                 rownames(distance_matrix)[j]][1])^2
                                     + (KEGG_mappings$Ycoord[KEGG_mappings$entryID ==
                                                                 colnames(distance_matrix)[k]][1]
                                        - KEGG_mappings$Ycoord[KEGG_mappings$entryID ==
                                                                   rownames(distance_matrix)[j]][1])^2)
                            # print(value)
                            distance_matrix[j,k] <- value
                        }
                    }
                    location <- which(distance_matrix == min(distance_matrix), arr.ind = TRUE)
                    entry_sources$entry1[i] <- colnames(distance_matrix)[location[1,2]]
                    entry_sources$entry2[i] <- rownames(distance_matrix)[location[1,1]]
                }
                if(length(entry1_locations) == 1 & length(entry2_locations) == 1){
                    entry_sources$entry1[i] <- entry1_locations
                    entry_sources$entry2[i] <- entry2_locations
                }
            }
        }
        entry_targets <- entry_edges[entry_edges$entry1accession == entry_accession &
                                         entry_edges$entry2accession %in% KEGG_mappings$entryACCESSION,]
        entry_targets <- entry_targets[!duplicated(entry_targets[,6:16]),]
        if(nrow(entry_targets) > 0){
            for(i in 1:nrow(entry_targets)){
                entry1_locations <- KEGG_mappings$entryID[KEGG_mappings$entrySYMBOL == entry_targets$entry1symbol[i]]
                entry2_locations <- KEGG_mappings$entryID[KEGG_mappings$entrySYMBOL == entry_targets$entry2symbol[i]]
                if(length(entry1_locations) > 1 | length(entry2_locations > 1)){
                    
                    distance_matrix <- matrix(nrow = length(entry2_locations),
                                              ncol = length(entry1_locations))
                    colnames(distance_matrix) <- entry1_locations
                    rownames(distance_matrix) <- entry2_locations
                    for(j in 1:nrow(distance_matrix)){
                        for (k in 1:ncol(distance_matrix)){
                            value <-
                                sqrt((KEGG_mappings$Xcoord[KEGG_mappings$entryID ==
                                                               colnames(distance_matrix)[k]][1]
                                      - KEGG_mappings$Xcoord[KEGG_mappings$entryID ==
                                                                 rownames(distance_matrix)[j]][1])^2
                                     + (KEGG_mappings$Ycoord[KEGG_mappings$entryID ==
                                                                 colnames(distance_matrix)[k]][1]
                                        - KEGG_mappings$Ycoord[KEGG_mappings$entryID ==
                                                                   rownames(distance_matrix)[j]][1])^2)
                            # print(value)
                            distance_matrix[j,k] <- value
                        }
                    }
                    location <- which(distance_matrix == min(distance_matrix), arr.ind = TRUE)
                    entry_targets$entry1[i] <- colnames(distance_matrix)[location[1,2]]
                    entry_targets$entry2[i] <- rownames(distance_matrix)[location[1,1]]
                }
                if(length(entry1_locations) == 1 & length(entry2_locations) == 1){
                    entry_targets$entry1[i] <- entry1_locations
                    entry_targets$entry2[i] <- entry2_locations
                }
            }
        }
        new_edges <- rbind(entry_sources, entry_targets)
        if(nrow(new_edges) > 0) {
            new_edges$primary_pathwayID <- strsplit(pathway_KGML@pathwayInfo@name, ":")[[1]][2]
            new_edges$primary_pathway <- pathway_KGML@pathwayInfo@title
            edges <- rbind(edges, new_edges)
        }
        print(p)
    }
    
    edges$edge_entry_ID <- paste0(edges$entry1, ":", edges$entry2)
    edge_ID_sets <- edges$edge_entry_ID[!duplicated(edges$edge_entry_ID)]
    edge_ID_mapper <- data.frame("edges" = edge_ID_sets, 
                                 edgeID = 1:length(edge_ID_sets),
                                 stringsAsFactors = FALSE)
    for(i in 1:nrow(edges)){
        edges$edgeID[i] <- edge_ID_mapper$edgeID[match(edges$edge_entry_ID[i], 
                                                       edge_ID_mapper$edges)]
    }
    edges <- edges[,names(edges) != "edge_entry_ID"]
    return(edges)
}



# KGML <- get_KGML("hsa04150")
# KEGG_mappings <- expand_KEGG_mappings(KGML)
# edges <- expand_KEGG_edges(KGML, KEGG_mappings)
# 

# 
# keggerize_edges <- function(entry_accession, KGML, KEGG_mappings, edges){
#     if(!"primary_pathway" %in% names(edges)){
#         edges$primary_pathwayID <- strsplit(KGML@pathwayInfo@name, ":")[[1]][2]
#         edges$primary_pathway <- KGML@pathwayInfo@title
#     }
#     entry_code <- paste0("hsa:",entry_accession)
#     entry_info <- keggGet(entry_code)
#     entry_pathways <- entry_info[[1]]["PATHWAY"]
#     entry_pathways <- names(entry_pathways$PATHWAY[1:length(unlist(entry_pathways))])
#     current_pathway <- strsplit(KGML@pathwayInfo@name, ":")[[1]][2]
#     entry_pathways <- entry_pathways[!entry_pathways %in% current_pathway]
#     for(i in 1:length(entry_pathways)){
#         pathway_KGML <- get_KGML(entry_pathways[i])
#         pathway_KEGG_mappings <- expand_KEGG_mappings(pathway_KGML)
#         pathway_edges <- expand_KEGG_edges(pathway_KGML, pathway_KEGG_mappings)
#         entry_edges <- pathway_edges[pathway_edges$entry1accession == entry_accession |
#                                          pathway_edges$entry2accession == entry_accession, ]
#         entry_sources <- entry_edges[entry_edges$entry1accession != entry_accession &
#                                          entry_edges$entry1accession %in% KEGG_mappings$entryACCESSION,]
#         # for(i in 1:nrow(entry_sources)){
#         #     entry1_locations <- KEGG_mappings$entryID[KEGG_mappings$entrySYMBOL == entry_sources$entry1symbol[i]]
#         #     entry2_locations <- KEGG_mappings$entryID[KEGG_mappings$entrySYMBOL == entry_sources$entry2symbol[i]]
#         #     distance_matrix <- matrix(nrow = length(entry2_locations), 
#         #                               ncol = length(entry1_locations))
#         #     colnames(distance_matrix) <- entry1_locations
#         #     rownames(distance_matrix) <- entry2_locations
#         #     for(j in 1:nrow(distance_matrix)){
#         #         for (k in 1:ncol(distance_matrix)){
#         #             value <- 
#         #                 sqrt((KEGG_mappings$Xcoord[KEGG_mappings$entryID ==
#         #                                               colnames(distance_matrix)[k]][1]
#         #                      - KEGG_mappings$Xcoord[KEGG_mappings$entryID ==
#         #                                                 rownames(distance_matrix)[j]][1])^2
#         #                      + (KEGG_mappings$Ycoord[KEGG_mappings$entryID ==
#         #                                                  colnames(distance_matrix)[k]][1]
#         #                         - KEGG_mappings$Ycoord[KEGG_mappings$entryID ==
#         #                                                    rownames(distance_matrix)[j]][1])^2)
#         #             print(value)
#         #             distance_matrix[j,k] <- value
#         #         }
#         #     }
#         #     
#         # 
#         # }
#         # for(i in 1:nrow(entry_sources)){
#         #     entry_sources$entry2[i] <- list(KEGG_mappings$entryID[KEGG_mappings$entrySYMBOL == entry_sources$entry2symbol[i]])
#         # }
#         # }
#         entry_targets <- entry_edges[entry_edges$entry1accession == entry_accession &
#                                          entry_edges$entry2accession %in% KEGG_mappings$entryACCESSION,]
#         new_edges <- rbind(entry_sources, entry_targets)
#         if(nrow(new_edges) > 0) {
#             new_edges$primary_pathwayID <- strsplit(pathway_KGML@pathwayInfo@name, ":")[[1]][2]
#             new_edges$primary_pathway <- pathway_KGML@pathwayInfo@title
#             
#             max_edge_ID <- max(edges$edgeID)
#             for(j in 1:nrow(new_edges)){
#                 new_edges$edgeID[j] <- max_edge_ID + j
#             }
#             edges <- rbind(edges, new_edges)
#         }
#         print(i)
#     }
#     edges$edge_entry_ID <- paste0(edges$entry1, ":", edges$entry2)
#     edge_ID_sets <- edges$edge_entry_ID[!duplicated(edges$edge_entry_ID)]
#     edge_ID_mapper <- data.frame("edges" = edge_ID_sets, 
#                                  edgeID = 1:length(edge_ID_sets),
#                                  stringsAsFactors = FALSE)
#     for(i in 1:nrow(edges)){
#         edges$edgeID[i] <- edge_ID_mapper$edgeID[match(edges$edge_entry_ID[i], 
#                                                        edge_ID_mapper$edges)]
#     }
#     edges <- edges[,names(edges) != "edge_entry_ID"]
#     return(edges)
# }
# 
# mtor_plus_mtor <- keggerize_edges(entry_accession = entry_accession, KGML = KGML,
#                                   KEGG_mappings = KEGG_mappings,
#                                   edges = edges)
# #dim(edges)
# #[1] 581  16
# # dim(mtor_plus_mtor)
# # [1] 757  18
# mtor_plus_mtor_plus_deptor <- keggerize_edges(entry_accession = "64798", KGML = KGML,
#                 KEGG_mappings = KEGG_mappings,
#                 edges = mtor_plus_mtor)
# 
# 
# 
# 
# 
# na_sources <- KEGG_gene_edges[is.na(KEGG_gene_edges$source),]
# na_targets<- KEGG_gene_edges[is.na(KEGG_gene_edges$target),]
# na_sources_and_targets <-  KEGG_gene_edges[is.na(KEGG_gene_edges$target)& is.na(KEGG_gene_edges$source),]
# 
# edges_complete_targets <- KEGG_gene_edges[!is.na(KEGG_gene_edges$target),] 
# complete_edges <- edges_complete_targets[!is.na(edges_complete_targets$source),]
# 
# # > dim(na_sources) + dim(na_targets) - dim(na_sources_and_targets)
# # [1] 2176   10
# # 
# # > dim(KEGG_gene_edges) - dim(complete_edges)
# # [1] 2176    0
# 
# 
# for(j in 1:nrow(distance_matrix)){
#     for (k in 1:ncol(distance_matrix)){
#         distance_matrix[j,k] <
#             sqrt((KEGG_mappings$Xcoord[KEGG_mappings$entryID ==
#                                            colnames(distance_matrix)[1]][1]
#                   - KEGG_mappings$Xcoord[KEGG_mappings$entryID ==
#                                              rownames(distance_matrix)[1]][1])^2
#                  + (KEGG_mappings$Ycoord[KEGG_mappings$entryID ==
#                                              colnames(distance_matrix)[1]][1]
#                     - KEGG_mappings$Ycoord[KEGG_mappings$entryID ==
#                                                rownames(distance_matrix)[1]][1])^2)
#     }
# }
# 
# KEGG_gene_edges <- complete_edges
# mtor_edges <- KEGG_gene_edges[KEGG_gene_edges$source == "MTOR" | KEGG_gene_edges$target == "MTOR",]
# 
# unique_mtor_edges <- mtor_edges[!duplicated(mtor_edges),]
# mtor_edges_all_pathways <- unique_mtor_edges
# 
# mtor_sources <- mtor_edges_all_pathways$source
# mtor_sources <- mtor_sources[!duplicated(mtor_sources)]
# # length(mtor_sources)
# # [1] 33
# mtor_targets <- mtor_edges_all_pathways$target
# mtor_targets <- mtor_targets[!duplicated(mtor_targets)]
# # length(mtor_targets)
# # [1] 26
# 
# mtor_KGML <- get_KGML("hsa04150")
# mtor_KEGG_mappings <- expand_KEGG_mappings(mtor_KGML)
# mtor_edges <- expand_KEGG_edges(mtor_KGML, mtor_KEGG_mappings)
# 
# sources_in_mtor_pathway <-  mtor_KEGG_mappings[mtor_KEGG_mappings$entrySYMBOL %in% mtor_sources,]
# # dim(sources_in_mtor_pathway)
# # [1] 28 16
# targets_in_mtor_pathway <-  mtor_KEGG_mappings[mtor_KEGG_mappings$entrySYMBOL %in% mtor_targets,]
# # dim(targets_in_mtor_pathway)
# # [1] 19 16
# internal_mtor_edges <- mtor_edges_all_pathways[mtor_edges_all_pathways$pathway_code == "path:hsa04150",]
