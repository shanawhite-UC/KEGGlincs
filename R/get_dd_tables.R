#' Import disease/drug tables from KEGG 
#' @description Get data tables for disease/drug information associated with 
#' selected pathway
#' @param pathwayid A KEGG pathway ID of the form "hsa12345" 
#' (only human pathways currently)
#' @return A data.frame object with either disease or drug information
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
