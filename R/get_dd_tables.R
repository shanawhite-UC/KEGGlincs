#' Import disease/drug tables from KEGG 
#' @description Get data tables for disease/drug information associated with 
#' selected pathway
#' @param pathwayid A KEGG pathway ID of the form "hsa12345" 
#' (only human pathways currently)
#' @return A data.frame object with either disease or drug information
#' @import XML
#' @export
#' @examples
#' RA_drug_table <- get_drug_table("hsa05323")

get_drug_table <- function(pathwayid){
    
    url <- paste0("http://www.kegg.jp/kegg-bin/pathway_dd_list?map=",
                  pathwayid)
    dd_table <- data.frame(
                XML::readHTMLTable(url, header = T, which = 4, as.data.frame = FALSE), 
                stringsAsFactors = FALSE)
    if (nrow(dd_table) > 0){
        d_table <- subset(dd_table, substring(dd_table$Disease.drug, 1, 1) == "D")
        if(nrow(d_table) == 0){
            warning("No associated drug targets in selected pathway")
            return()
        }
        names(d_table) <- c("drug_KEGG_ID", "drug_name", "gene_target")
        d_table$gene_target <- strsplit(d_table$gene_target, " ")

        for(i in 1:nrow(d_table)){
            l <- length(unlist(d_table$gene_target[i]))
            d_table$drug_KEGG_ID[i] <- list(rep(d_table$drug_KEGG_ID[i], l))
            d_table$drug_name[i] <- list(rep(d_table$drug_name[i], l))
        }
        
        long_drug <- data.frame("drug_KEGG_ID" = unlist(d_table$drug_KEGG_ID), 
                                "drug_name" = unlist(d_table$drug_name), 
                                "gene_target" = unlist(d_table$gene_target),
                                 stringsAsFactors = FALSE)
        for (i in 1:nrow(long_drug)){
            long_drug$gene_id[i] <- strsplit(long_drug$gene_target[i], "\\(")[[1]][1]
            long_drug$gene_symbol[i] <- regmatches(long_drug$gene_target[i], gregexpr("(?<=\\().*?(?=\\))", long_drug$gene_target[i], perl=T))[[1]]
        }
        drops <- "gene_target"
        d_table <- long_drug[, names(long_drug) != drops]
        
        return(d_table)
    }
    else {
        warning("No drugs or diseases associated with selected pathway")
        return()
    }
}
#' @rdname get_drug_table 
get_disease_table <- function(pathwayid){
    
    url <- paste0("http://www.kegg.jp/kegg-bin/pathway_dd_list?map=",
                  pathwayid)
    dd_table <- data.frame(
        XML::readHTMLTable(url, header = T, which = 4, as.data.frame = FALSE), 
        stringsAsFactors = FALSE)
    if (nrow(dd_table) > 0){
        d_table <- subset(dd_table, substring(dd_table$Disease.drug, 1, 1) == "H")
        if(nrow(d_table) == 0){
            warning("No diseases associated with selected pathway")
            return()
        }
        return(d_table)
    }
    else {
        warning("No drugs or diseases associated with selected pathway")
        return()
    }
}

#' @examples
#' RA_disease_table <- get_disease_table("hsa05323")
