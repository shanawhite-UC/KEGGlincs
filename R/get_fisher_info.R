#' Perform Fisher's Exact test for edges in pathway
#' @description Obtain a measure for strength and significance for the 
#' relationship (i.e. an edge) based on the concordance/discordance of 
#' UP-and-DOWN regulated genes shared by two different experimental 
#' gene-knockouts
#' Intended for use within \code{\link{overlap_info}}
#' @export
#' @importFrom stats fisher.test p.adjust p.adjust.methods
#' @param edges The set of eges to be analyzed; Although the intended use is for
#' LINCS data overlaps, the function should work with any typical data object as
#' long as it has columns labeled ("UP", "DOWN", "UK1_DK2", "DK1_UK2") that 
#' contain integer values.
#' @param method The method to correct/adjust p-values for multiple testing. 
#' For available methods, type 'p.adjust.methods' into command promt and 
#' press enter.
#' @return The input edge data.frame object with additional columns containing
#' the results of the applied statistical test  
#' @examples
#' ex.data <- data.frame("UP" = c(70,6), "DOWN" = c(8,20),
#'                     "UK1_DK2" = c(4,47), "DK1_UK2" = c(3, 28))
#' 
#' overlaps <- get_fisher_info(ex.data, method = "BH")

get_fisher_info <-
function(edges, method){
    if (!method %in% p.adjust.methods){
        warning("Method for p-value adjustment is not valid. For valid 
                adjustment methods please type 'p.adjust.methods' in the 
                command line and press enter. P-values will not be adjusted")
        method = "none"
    }
    for (i in 1:nrow(edges)){
        fisher_test<- fisher.test(matrix(c(edges$UP[i] + 1,edges$UK1_DK2[i] + 1,
                                edges$DK1_UK2[i] + 1 ,edges$DOWN[i] + 1), 2,2))
        edges$fisher_OR[i] <- unlist(fisher_test[[3]])
        edges$fisher_pval[i] <- unlist(fisher_test[[1]])
        if (edges$fisher_OR[i] <= 1) {
            edges$fisher_sign[i] <- -1
        }
        else {
            edges$fisher_sign[i] <- 1
        }
        edges$summary_score[i] <- 
            -log(edges$fisher_pval[i])*edges$fisher_sign[i]
    }
    raw_pvals <- edges$fisher_pval
    adj_pvals <- p.adjust(raw_pvals, method = method)
    edges$p_value <-  adj_pvals
    for (i in 1:nrow(edges)){
        if (edges$p_value[i] <= 0.05){
            edges$significant[i] <- 1
        }
        else {
            edges$significant[i] <- 0
        }
    }
    return(edges)
}
