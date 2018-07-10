#' Quickly access LINCS L1000 CGS's for set of perturbagens (KO's)
#' @description Translate raw CGS data to easy-to-use format
#' @export
#' @param cell_line Choose from the set of cell lines: 
#' (A375,A549,ASC,HA1E,HCC515,HEK293T,HEKTE,HEPG2,HT29,MCF7,NCIH716,NPC,PC3,
#' SHSY5Y,SKL,SW480,VCAP)  
#' @param edge_id The numeric value for the edge_id
#' @param data_type Choose from data types: (100_full, 100_bing, 50_lm)
#' @param pert_time Choose from (6,24,48,96,120,144,168)
#' @param pathway_nodes Keep NA unless certain set of perturbagens is designated
#' @return A data frame with conveniently formatted LINCS L100 CGS
#' @examples 
#' MCF_LM_50 <- grab_KO_data("MCF7")
#' MCF_BING_100 <- grab_KO_data("MCF7", data_type = "100_bing")
#' 
grab_KO_data <- function(cell_line, pert_time = 96, data_type = "50_lm", 
                         pathway_nodes = NA){
    data("KO_data", envir = environment())
    KO_data <- get("KO_data")
    keeps <- c(names(KO_data[1:3]), paste0("up",data_type), 
               paste0("dn",data_type))    
    suppressWarnings(    
        if(is.na(pathway_nodes)){
            KO_by_CT <- KO_data[KO_data$pert_time == pert_time & 
                                    KO_data$cell_id == cell_line, keeps]
        }
    )
    suppressWarnings(
        if(!is.na(pathway_nodes)){
            KO_by_CT <- KO_data[KO_data$pert_time == pert_time & 
                                    KO_data$cell_id == cell_line &
                                    KO_data$pert_desc %in% pathway_nodes, keeps]
        }
    )
    names(KO_by_CT)[c(4,5)] <- c("up", "down")
    KO_by_CT <- KO_by_CT[,c(2,4:5)]
    data("conversion_key")
    data("L1000_LM_genes")
    
    for (i in 1:nrow(KO_by_CT)){
        
        KO_by_CT$up_SYMBOL[i] <- 
            list(conversion_key$pr_gene_symbol[
                which(conversion_key$pr_id %in% 
                unlist(strsplit(KO_by_CT$up[i], ";")))])
        KO_by_CT$down_SYMBOL[i] <-
            list(conversion_key$pr_gene_symbol[
                which(conversion_key$pr_id %in% 
                unlist(strsplit(KO_by_CT$down[i], ";")))])
        KO_by_CT$up_count[i] <- length(unlist(KO_by_CT$up_SYMBOL[i]))
        KO_by_CT$down_count[i] <- length(unlist(KO_by_CT$down_SYMBOL[i]))
    }
    KO_by_CT <- KO_by_CT[,c(1,4,5)]
    names(KO_by_CT) <- c("knockout", "up", "down")
    rownames(KO_by_CT) <- KO_by_CT$knockout
    return(KO_by_CT)
}