#' Add edges from disease/drug tables  
#' @description Expand edge mappings to include drugs/drug targets
#' for selected pathway
#' @param edges A data.frame object obtained by using the function `expand_kegg_edges`
#' @param KEGG_mappings A data.frame object obtained by using the function `expand_kegg_mappings`
#' @param kegg_drug_table A data.frame object obtained by using the function `get_drug_table`
#' @return A data.frame object similar to the expanded edges data frame but with additional
#' edges representing known drugs/drug targets
#' @export
#' @examples

#' end_res_KGML <- get_KGML("hsa01522")
#' end_res_KEGG_mappings <- expand_KEGG_mappings(end_res_KGML)
#' end_res_edges <- expand_KEGG_edges(end_res_KGML, end_res_KEGG_mappings)

#' end_res_drugs <- get_drug_table("hsa01522")

#' edges_plus_kdrug <- add_kegg_drugs(end_res_edges, end_res_KEGG_mappings, end_res_drugs)

add_kegg_drugs <- function(edges, KEGG_mappings, kegg_drug_table){

    edges$u_ID <- paste0(edges$entry1accession,":", edges$entry2accession)
    kegg_drug_table$u_ID <- paste0(kegg_drug_table$drug_KEGG_ID,":", kegg_drug_table$gene_id)
    
    drugs_to_add <- subset(kegg_drug_table, !kegg_drug_table$u_ID %in% edges$u_ID)
    edges_to_add <- data.frame(
       "edgeID" = seq(from = (1 + max(edges$edgeID)), to = (nrow(drugs_to_add) + max(edges$edgeID))), 
       "entry1accession" = drugs_to_add$drug_KEGG_ID,
       "entry2accession" = drugs_to_add$gene_id,
       "entry1" = NA,
       "entry2" = NA,
       "type" = "PCrel",
       "subtype1" = "from_kegg_drug",
       "value" = "--",
       "subtype2" = NA,
       "value2" = NA,
       "specific_subtype" = "from_kegg_drug",
       "is_direct" = 1,
       "entry1type" = "compound",
       "entry2type" = "gene",
       "entry1symbol" = drugs_to_add$drug_name,
       "entry2symbol" = drugs_to_add$gene_symbol,
       "u_ID" = drugs_to_add$u_ID,
       stringsAsFactors = FALSE
    )
    for (i in 1:nrow(edges_to_add)){
        edges_to_add$entry2[i] <- KEGG_mappings$entryID[KEGG_mappings$entryACCESSION == 
                                                            edges_to_add$entry2accession[i]][1]
        if (edges_to_add$entry1accession[i] %in% KEGG_mappings$entryACCESSION){
            edges_to_add$entry1[i] <- KEGG_mappings$entryID[KEGG_mappings$entryACCESSION == 
                                                                edges_to_add$entry1accession[i]][1]
        }
    }
    # Decide whether or not to map to all nodes of gene on map
    # for (i in 1:nrow(edges_to_add)){
    #     edges_to_add$entry2[i] <- list(KEGG_mappings$entryID[KEGG_mappings$entryACCESSION == edges_to_add$entry2accession[i]])
    # }
    all_edges <- rbind(edges[, -c(17)], edges_to_add[, -c(17)])

}

#load("/opt/raid10/genomics/shana/KEGG_drug_targets/via_cts_convert/compound_cgs_KEGG_filtered.rda")
