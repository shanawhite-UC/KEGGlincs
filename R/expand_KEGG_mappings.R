#' Get detailed KEGG mapping information for each map entity
#' @description Extract mapping information from KGML object and normalize 
#' mappings based on multi-valued name attribute
#'
#' @param KGML_file An object of formal class KEGGPathway
#' @param convert_KEGG_IDs A logical indicator; if set to FALSE will run faster
#'  however genes and compounds will remain labeled via KEGG codes (compounds) 
#'  or accession numbers (genes).  This option must be taken into account if 
#'  data is being added.  For example, the genes in 'KO_data' are identified by
#'  symbols, thus it is neccessary to retain the default option to convert IDs
#'  to symbols when planning to add edge data of this type.
#'
#' @return A dataframe object with unique entry information for all [node]
#'  objects documented in the KEGG pathway. 
#'         Note that if mutiple objects (i.e. genes or compounds) have the same 
#'         entryID, this indicates that they share the same node [location] 
#'         in the pathway.  
#' @export
#' @import org.Hs.eg.db 
#' @import KEGGREST
#' @import AnnotationDbi
#' @importMethodsFrom KEGGgraph edges
#' @importMethodsFrom KEGGgraph nodes
#' @examples
#' p53_KGML <- get_KGML("hsa04115")
#' p53_KEGG_mappings <-  expand_KEGG_mappings(p53_KGML, FALSE)

expand_KEGG_mappings <-
function (KGML_file, convert_KEGG_IDs = TRUE){
    nodesN <- KEGGgraph::nodes(KGML_file)
    entryMAP<- data.frame(entry = seq(1:length(nodesN)))
    for (i in 1:nrow(entryMAP)){
        entryMAP$entryID[i] <- nodesN[[i]]@entryID
        entryMAP$entryTYPE[i] <- nodesN[[i]]@type
        entryMAP$entryIDENTIFIER[i] <- list(nodesN[[i]]@name)
        if (entryMAP$entryTYPE[i] == "group"){
            entryMAP$groupMEMBERS[i] <-
                list(nodesN[[entryMAP$entryID[i]]]@component)
            entryMAP$entryIDENTIFIER[i] <- 
                list(entryMAP$entryIDENTIFIER[which(entryMAP$entryID %in% 
                                            unlist(entryMAP$groupMEMBERS[i]))])
        }
        else {
            entryMAP$groupMEMBERS[i] <- NA
        }
        if (entryMAP$entryTYPE[i] == "gene" | 
            entryMAP$entryTYPE[i] == "compound"|
            entryMAP$entryTYPE[i] == "group"){
                entryMAP$entryACCESSION[i]<- list(gsub(".*:", "", 
                                        unlist(entryMAP$entryIDENTIFIER[i])))
        }
        else {
            entryMAP$entryACCESSION[i] <- NA
        }
        entryMAP$entryTYPE[i] <- nodesN[[i]]@type
        entryMAP$entryNAMES[i] <- nodesN[[i]]@graphics@name
        entryMAP$FGcolor[i] <- nodesN[[i]]@graphics@fgcolor
        entryMAP$BGcolor[i] <- nodesN[[i]]@graphics@bgcolor
        entryMAP$shape[i] <- nodesN[[i]]@graphics@type
        entryMAP$Xcoord[i] <- nodesN[[i]]@graphics@x
        entryMAP$Ycoord[i] <- nodesN[[i]]@graphics@y
        entryMAP$width[i] <- nodesN[[i]]@graphics@width
        entryMAP$height[i] <- nodesN[[i]]@graphics@height
        if(is.na(entryMAP$entryNAMES[i])){
            entryMAP$entryNAMES[i] <- "NoNameInKGML"
        }
    }
    entryMAP <- entryMAP[,-c(1,4)]
    entryMAP_groupings <- entryMAP[entryMAP$entryTYPE == "group",]
    entryMAP_extra <- entryMAP[entryMAP$entryTYPE == "map",]
    entryMAP_extra$entrySYMBOL <- NA
    entryMAP_dictionary <- entryMAP[entryMAP$entryTYPE == "gene" | 
                                    entryMAP$entryTYPE == "compound",]
    for (i in 1:nrow(entryMAP_dictionary)){  
        if (entryMAP_dictionary$entryTYPE[i] == "compound"){
            entryMAP_dictionary$entryACCESSION[i] <- 
            entryMAP_dictionary$entryACCESSION[[i]][1]
        }
    }
    if(convert_KEGG_IDs){
        for (i in 1:nrow(entryMAP_dictionary)){  
            if (entryMAP_dictionary$entryTYPE[i] == "compound"){
                testval <- 
                    keggGet(entryMAP_dictionary$entryACCESSION[i])[[1]]$NAME
                if(!is.null(testval)){
                    entryMAP_dictionary$entryNAMES[i] <- list(testval)
                    entryMAP_dictionary$entryNAMES[i] <- 
                        list(gsub(";", "",
                                unlist(entryMAP_dictionary$entryNAMES[i])))
                }
            }
        }
    }
    for (i in 1:nrow(entryMAP_dictionary)){  
        l = length(unlist(entryMAP_dictionary$entryACCESSION[i]))
        entryMAP_dictionary$entryID[i] <- 
            list(rep(entryMAP_dictionary$entryID[i], l))
        entryMAP_dictionary$entryTYPE[i] <- 
            list(rep(entryMAP_dictionary$entryTYPE[i], l))
    }
    emd_temp <- data.frame(cbind(unlist(entryMAP_dictionary$entryID), 
                                unlist(entryMAP_dictionary$entryTYPE), 
                                unlist(entryMAP_dictionary$entryACCESSION)), 
                                stringsAsFactors = FALSE)
    names(emd_temp) <- names(entryMAP[c(1,2,4)])
    for (i in 1:nrow(emd_temp)){
        if (emd_temp$entryTYPE[i] == "gene"){
            if (convert_KEGG_IDs){
                emd_temp$entrySYMBOL[i] <- 
                suppressMessages(as.character(AnnotationDbi::select(
                org.Hs.eg.db::org.Hs.eg.db,keys = emd_temp$entryACCESSION[i], 
                columns=c("SYMBOL"),keytype="ENTREZID")[2]))
            }
            else {
                emd_temp$entrySYMBOL[i] <- emd_temp$entryACCESSION[i]
            }
        }
        else if (emd_temp$entryTYPE[i] == "compound"){
            emd_temp$entrySYMBOL[i] <- unlist(
                entryMAP_dictionary$entryNAMES[
                    entryMAP_dictionary$entryACCESSION == 
                    emd_temp$entryACCESSION[i]])[1]
        }
        else {
            emd_temp$entrySYMBOL[i] <- NA
        }
    }
    for (i in 1:nrow(emd_temp)){
        if (gsub("[[:digit:]]","",emd_temp$entryACCESSION[i]) == "G" & 
            is.na(emd_temp$entrySYMBOL[i])){
            emd_temp$entrySYMBOL[i] <- emd_temp$entryACCESSION[i]
        }
    }
    emd_temp$groupMEMBERS <- NA
    for(i in 1:nrow(entryMAP_dictionary)){
        entryMAP_dictionary$entryID[i] <- 
            unlist(entryMAP_dictionary$entryID[i])[1] 
        entryMAP_dictionary$entryTYPE[i] <- 
            unlist(entryMAP_dictionary$entryTYPE[i])[1]
        if(entryMAP_dictionary$entryTYPE[i] == "compound"){
            entryMAP_dictionary$entryNAMES[i] <- 
                toString(unlist(entryMAP_dictionary$entryNAMES[[i]]))
        }
    }
    entryMAP_dictionary <- entryMAP_dictionary[,-c(2:4)]
    entryMAP_dictionary <- merge(emd_temp, entryMAP_dictionary, by = "entryID")
    if (nrow(entryMAP_groupings) > 0){
        for(i in 1:nrow(entryMAP_groupings)){
            entryMAP_groupings$entrySYMBOL[i] <- 
                list(entryMAP_dictionary$entrySYMBOL[match(unlist(
                    entryMAP_groupings$entryACCESSION[i]),
                    entryMAP_dictionary$entryACCESSION)])
        }
        entryMAP_dictionary <- 
            rbind(entryMAP_extra, entryMAP_dictionary, entryMAP_groupings)
        no_groups <- FALSE
    }
    if (nrow(entryMAP_groupings) == 0){
        entryMAP_dictionary <- rbind(entryMAP_extra, entryMAP_dictionary)
        entryMAP_dictionary$groupMEMBERS <- NA
        entryMAP_dictionary$groupID <- NA
        no_groups = TRUE
    }

    entryMAP_dictionary <- 
        entryMAP_dictionary[order(as.double(entryMAP_dictionary$entryID)),]
    rownames(entryMAP_dictionary) <- seq(1:nrow(entryMAP_dictionary))
    entryMAP_dictionary$entryNAMES <- gsub("...", "", 
                                            entryMAP_dictionary$entryNAMES,
                                            fixed = TRUE)
    for (i in 1:nrow(entryMAP_dictionary)){
        if (entryMAP_dictionary$entryTYPE[i] =="gene"){
            entryMAP_dictionary$LABEL[i] <- strsplit(as.character(
                entryMAP_dictionary$entryNAMES[[i]]), ",")[[1]][1]
        }
        else if (entryMAP_dictionary$entryTYPE[i] == "compound"){
            entryMAP_dictionary$LABEL[i] <- entryMAP_dictionary$entrySYMBOL[i]
        }
        else {
            entryMAP_dictionary$LABEL[i] = entryMAP_dictionary$entryNAMES[i]
        }
    }
    for (i in 1:nrow(entryMAP_dictionary)){
        if (entryMAP_dictionary$entryTYPE[i] == "group"){
            entryMAP_dictionary$entryNAMES[i] <- 
                paste0(paste(unique(entryMAP_dictionary$LABEL[which(
                entryMAP_dictionary$entryID %in% 
                unlist(entryMAP_dictionary$groupMEMBERS[i]))]),collapse = ":"),
                " Complex")
            entryMAP_dictionary$LABEL[i] <- entryMAP_dictionary$entryNAMES[i]
        }
    }
    if (!no_groups){
        for (i in 1:nrow(entryMAP_groupings)){
            l <- length(unlist(entryMAP_groupings$groupMEMBERS[i]))
            entryMAP_groupings$entryID[i] <-  
                list(rep(entryMAP_groupings$entryID[i], l))
        }
        group_mapper <- data.frame("groupID" = 
                                    unlist(entryMAP_groupings$entryID), 
                            "entryID" = 
                                unlist(entryMAP_groupings$groupMEMBERS),
                            stringsAsFactors = FALSE)
        for (i in 1:nrow(entryMAP_dictionary)){
            if (entryMAP_dictionary$entryID[i] %in% group_mapper$entryID){
                entryMAP_dictionary$groupID[i] <-
                    group_mapper$groupID[match(entryMAP_dictionary$entryID[i],
                                        group_mapper$entryID)]
            }
            else {
                entryMAP_dictionary$groupID[i] <- NA
            }
        }
    }

    num_edges <- length(KGML_file@edges)
    if (num_edges > 0){
        map_edge_data<- KEGGgraph::edges(KGML_file)
        edges<- data.frame(edgeID = seq(1:length(map_edge_data )))
        for (i in 1:nrow(edges)){
            edges$entry1[i] <- as.character(map_edge_data[[i]]@entry1ID)
            edges$entry2[i] <- as.character(map_edge_data[[i]]@entry2ID)
        }
        for (i in 1:nrow(entryMAP_dictionary)){
            if (entryMAP_dictionary$entryID[i] %in% edges$entry1 | 
                entryMAP_dictionary$entryID[i] %in% edges$entry2){
                    entryMAP_dictionary$in_relationship[i] = 1
            }
            else if (entryMAP_dictionary$groupID[i] %in% edges$entry1 | 
                    entryMAP_dictionary$groupID[i] %in% edges$entry2){
                    entryMAP_dictionary$in_relationship[i] = 1
            }
            else {
                entryMAP_dictionary$in_relationship[i] = 0
            }
        }
    }
    if (num_edges == 0){
      entryMAP_dictionary$in_relationship = 0
    }
    for (i in 1:nrow(entryMAP_dictionary)){
        if (is.na(entryMAP_dictionary$entrySYMBOL[i])){
            entryMAP_dictionary$entrySYMBOL[i] <- entryMAP_dictionary$LABEL[i]
        }
    }
    map_nodes <- entryMAP_dictionary$LABEL[entryMAP_dictionary$entryTYPE=="map"]
    j = 0
    for (i in 1:length(map_nodes)){
        if (substring(map_nodes[i],1,5) == "TITLE"){
            j = j+1
        }
    }
    if (j == 0){
        title <- KGML_file@pathwayInfo@title
        min_Y <- min(entryMAP_dictionary$Ycoord)
        min_X <- min(entryMAP_dictionary$Xcoord)+100
        title_node <- data.frame("entryID" = 0, "entryTYPE" = "map", 
                                "groupMEMBERS" = NA, "entryACCESSION" = NA, 
                                "entryNAMES" = title,
                                "FGcolor" = "#000000", "BGcolor" = "#FFFFFF",
                                "shape" = "rectangle", "Xcoord" = min_X, 
                                "Ycoord"= min_Y, "height" = 25, "width" = 220, 
                                "entrySYMBOL"= title, "groupID"=NA,
                                "LABEL"=paste0("TITLE:",title), 
                                "in_relationship" = 0)
        entryMAP_dictionary <- rbind(entryMAP_dictionary, title_node)
    }
    return(entryMAP_dictionary)
}