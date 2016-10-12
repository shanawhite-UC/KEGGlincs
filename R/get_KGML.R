#' Download and parse KGML file
#' @export  
#' @param pathwayid A KEGG pathway ID of the form "hsa12345" 
#' (only human pathways currently)
#' @param get_if_no_edges A logical indicator; if pathway has no edges 
#' returns null value if set to TRUE
#' @return an object of Formal class KEGGPathway
#' @examples 
#' mtor_KGML <- get_KGML("hsa04150")
#' 
#' # Some pathways contain only node information; since the purpose of this
#' # package is to explore pathways in an edge-focused manner, the default
#' # options return a warning message instead of a parsed KGML file if the 
#' # input pathway has no edges.
#' ribosome_KGML <- get_KGML("hsa03020") 
#' ribosome_KGML <- get_KGML("hsa03020", get_if_no_edges = TRUE) 
#' 

get_KGML <- 
function(pathwayid, get_if_no_edges = FALSE){
    tmp <-  suppressMessages(KEGGREST::keggGet(pathwayid, "kgml"))
    if(nchar(tmp) == 0){
        warning("Selected pathway does not have an associated KGML file")
        KGML_file <- NA
        return(KGML_file)
    }
    if(length(grep("relation",tmp)) == 0){
        if(get_if_no_edges){
            KGML_file <- KEGGgraph::parseKGML(tmp)
            return(KGML_file)
        }
        else {
            warning("Selected pathway does not contain any documented edges")
            return(NA)
        }
    }
    KGML_file <- KEGGgraph::parseKGML(tmp)
    return(KGML_file)
}
