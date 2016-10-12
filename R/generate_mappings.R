#' The 'boilerplate' for this package's desired graph style 
#' @description Generates an object that can be converted to a JSON file 
#' and subsequently applied to the graph for the markup specified by this
#'  package and the layout mirroring KEGG.
#' Intended for use within \code{\link{cyto_vis}}
#' @export
#' @param style_name An argument to name style; when used inside 
#' of \code{\link{cyto_vis}} no name is needed
#' @param map_edge_width A logical indicator; if FALSE no continuous mapping 
#' of edge width will be applied
#' @param edge_width_attribute The attribute that will be used for edge width; 
#' if data is not added or the attribute is not part of the graphing 
#' information, the edge width will default to 1. 
#' @param min_score The minimum attribute value for the column used to map 
#' edge width
#' @param max_score The maximum attribute value for the column used to map 
#' edge width
#' @return A list that can be converted to a JSON file to apply desired 
#' style/layout in Cytoscape  
#' @examples 
#' style.name = "myKEGGstyle"
#' mappings <- generate_mappings(style.name, FALSE)

generate_mappings <-
function(style_name, map_edge_width, edge_width_attribute, 
            min_score, max_score){
    style.name = style_name
    
    # Visual Mappings
    #######
    ##defaults
    def.node.font <- list(visualProperty = "NODE_LABEL_FONT_FACE", 
                            value = "TeXGyrePagella-Regular")
    
    defaults <- list(def.node.font)
    #######
    
    ####transparency
    pair1a = list(key = "group",value = "0")
    pair1b = list(key = "gene", value = "150")
    pair1c = list(key = "compound", value = "255")
    pair1d = list(key = "map",value = "255")
    discrete.mappings.1 = list(pair1a, pair1b, pair1c, pair1d)
    ####node shape
    pair2a = list(key = "circle", value = "Ellipse")
    pair2b = list(key = "roundrectangle", value = "Round_Rectangle")
    pair2c = list(key = "rectangle",  value = "Rectangle")
    discrete.mappings.2 = list(pair2a, pair2b, pair2c)
    ###edge target arrow shape
    pair3a = list(key = "activation", value = "Delta")
    pair3b = list(key = "expression", value = "Arrow")
    pair3c = list(key = "indirect_effect", value = "Diamond")
    pair3d = list(key = "inhibition",  value = "T")
    pair3e = list(key = "inhibition_phosphorylation",  value = "T")
    pair3f = list(key = "irreversible",  value = "Arrow")
    pair3g = list(key = "phosphorylation",  value = "Diamond")
    pair3h = list(key = "repression",  value = "T")
    discrete.mappings.3 = list(pair3a, pair3b, pair3c, pair3d, pair3e, pair3f, 
                                pair3g, pair3h)
###node label location
    pair4a = list(key = "compound", value = "S,N,c,0.0,0.0")
    pair4b = list(key = "group",  value = "E,W,c,2.0,0.0")
    discrete.mappings.4 = list(pair4a, pair4b)
###dashed edge for indirect effects
    pair5 = list(key = "0.0",  value = "Equal_Dash")
    discrete.mappings.5 = list(pair5)
###transparancy for edges [based on existence in data : 
###    not needed for 'default' graphs]
    pair6a = list(key = "0.0",value = "100")
    pair6b = list(key = "1.0", value = "255")
    discrete.mappings.6 = list(pair6a, pair6b)
    
    node.transparency = list(mappingType="discrete",
                                mappingColumn="entryTYPE",
                                mappingColumnType="String",
                                visualProperty="NODE_TRANSPARENCY",
                                map = discrete.mappings.1)
    node.shape = list(mappingType="discrete", 
                         mappingColumn="shape",
                         mappingColumnType="String",
                          visualProperty="NODE_SHAPE",
                         map = discrete.mappings.2)
    node.label.position = list(mappingType="discrete",
                               mappingColumn="entryTYPE",
                               mappingColumnType="String",
                               visualProperty="NODE_LABEL_POSITION",
                               map = discrete.mappings.4)
    node.color = list(mappingType="passthrough",
                        mappingColumn="BGcolor",
                        mappingColumnType="String",
                        visualProperty="NODE_FILL_COLOR")
    node.label = list(mappingType="passthrough",
                        mappingColumn="LABEL",
                        mappingColumnType="String",
                        visualProperty="NODE_LABEL")
    node.label.width = list(mappingType="passthrough",
                        mappingColumn="width",
                        mappingColumnType="String",
                        visualProperty="NODE_LABEL_WIDTH")
    node.Xcoord = list(mappingType="passthrough",
                       mappingColumn="Xcoord",
                       mappingColumnType="Double",
                       visualProperty="NODE_X_LOCATION")
    node.Ycoord = list(mappingType="passthrough",
                        mappingColumn="Ycoord",
                        mappingColumnType="Double",
                        visualProperty="NODE_Y_LOCATION")
    node.width = list(mappingType="passthrough",
                        mappingColumn="width",
                        mappingColumnType="Double",
                        visualProperty="NODE_WIDTH")
    node.height = list(mappingType="passthrough",
                        mappingColumn="height",
                        mappingColumnType="Double",
                        visualProperty="NODE_HEIGHT")
    node.border.width = list(mappingType="passthrough",
                                mappingColumn="border_width",
                                mappingColumnType="Double",
                                visualProperty="NODE_BORDER_WIDTH")
    node.label.font.size = list(mappingType="passthrough",
                                mappingColumn="label_font_size",
                                mappingColumnType="Double",
                                visualProperty="NODE_LABEL_FONT_SIZE")
    node.tooltip = list(mappingType="passthrough",
                        mappingColumn="entryNAMES",
                        mappingColumnType="String",
                        visualProperty="NODE_TOOLTIP")
    node.mappings = list(node.color, node.label, node.label.width, 
                        node.label.position, node.border.width,node.Xcoord, 
                        node.Ycoord, node.width, node.height, 
                        node.transparency, node.shape, node.label.font.size, 
                        node.tooltip)
    edge.target.arrow.shape = list(mappingType="discrete",
                                    mappingColumn="subtype1", 
                                    mappingColumnType="String",
                                    visualProperty="EDGE_TARGET_ARROW_SHAPE",
                                    map = discrete.mappings.3)
    edge.line.type = list(mappingType="discrete",
                            mappingColumn="is_direct",
                            mappingColumnType="Double",
                            visualProperty="EDGE_LINE_TYPE",
                            map = discrete.mappings.5)
    edge.transparency = list(mappingType="discrete",
                                mappingColumn="has_data",
                                mappingColumnType="Double",
                                visualProperty="EDGE_TRANSPARENCY",
                                map = discrete.mappings.6)
    edge.target.arrow.color = list(mappingType="passthrough",
                                    mappingColumn="color",
                                    mappingColumnType="String",
                                    visualProperty=
                                    "EDGE_TARGET_ARROW_UNSELECTED_PAINT")
    edge.label = list(mappingType="passthrough",
                        mappingColumn="edge_label",
                        mappingColumnType="String",
                        visualProperty="EDGE_LABEL")
    edge.color = list(mappingType="passthrough",
                        mappingColumn="color",
                        mappingColumnType="String",
                        visualProperty="EDGE_STROKE_UNSELECTED_PAINT")
    edge.tooltip = list(mappingType="passthrough",
                        mappingColumn="tooltip",
                        mappingColumnType="String",
                        visualProperty="EDGE_TOOLTIP")
    
    if (map_edge_width){
        point1 = list(value=min_score, lesser= "1", equal="1", greater="1")
        point2 = list(value=max_score, lesser="20", equal="20", greater="20")
        edge.width.continuous.points = list(point1, point2)
        edge.width = list(mappingType="continuous",
                            mappingColumn = edge_width_attribute,
                            mappingColumnType="Double",
                            visualProperty="EDGE_WIDTH",
                            points = edge.width.continuous.points)
        edge.mappings = list(edge.target.arrow.color, edge.target.arrow.shape, 
                            edge.line.type, edge.color, edge.tooltip, 
                            edge.label, edge.transparency, edge.width)
    }
    else {
        edge.mappings = list(edge.target.arrow.color, edge.target.arrow.shape, 
                            edge.line.type, edge.color,
                            edge.tooltip, edge.label)
    }
    mappings <- c(node.mappings, edge.mappings)
    return(list(defaults,mappings))
}

