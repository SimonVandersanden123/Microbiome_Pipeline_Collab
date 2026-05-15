#' @title forceGraph
#'
#' @description node-link visualization of the microbiome composition graph
#'
#' @param graph_object microbiome composition graph build with build_mc_graph().
#' @param source_color color variable to set color on source nodes (samples).
#' @param target_color color variable to set color on target nodes (ASVs).
#' @param groupingVariable grouping variable between the samples.
#' @param forceOnGroup add repulsion force on groups. Boolean.
#' @param nodeStroke color of node stroke. defualt is "black". custum mapping can be done by providing an array in the following format c('variable', 'scaleType', 'domain', 'range', 'unknown').
#' @param nodeStrokeWidth width of the node stroke. default is 1.5.
#' @param nodeRadius_source radius of source nodes. default is 10. custum mapping can be done by providing an array in the following format c('variable', 'scaleType', 'domain', 'range', 'unknown').
#' @param nodeRadius_target radius of target nodes. default is 4. custum mapping can be done by providing an array in the following format c('variable', 'scaleType', 'domain', 'range', 'unknown').
#' @param linkStrokeWidth width of lines between connected nodes. default is 1.5. custum mapping can be done by providing an array in the following format c('variable', 'scaleType', 'domain', 'range', 'unknown').
#' @param linkStroke color of lines between nodes. default is "gray". custum mapping can be done by providing an array in the following format c('variable', 'scaleType', 'domain', 'range', 'unknown').
#' @param linkStrokeOpacity opacity of lines between nodes. default is 0.6. custum mapping can be done by providing an array in the following format c('variable', 'scaleType', 'domain', 'range', 'unknown').
#' @param nodeOpacity opacity of nodes. default is 1.
#' @param nodeStrokeOpacity opacity of nodes' stroke. default is 1.
#' @param nodeTitle array of variables you want to display in the tooltip. c('variable1','variable2')
#' @param highlightAbundance show taxa abundance within sample when hovering over a sample node and show taxon abundance in samples when hovering over a taxon node. Boolean.
#' @param width width of the graph. default is full available width.
#' @param height hieght of the graph. default is full available height.
#'
#' @return interactive figure
#'
#' @examples
#' # Default
#' forceGraph(graph_object = mcg)
#' # Color source and target nodes
#' forceGraph(graph_object = mcg, source_color = "alphaShannon", target_color = "abundance")
#' # Set the node attributes to be printed while hovering on a node
#' forceGraph(graph_object = mcg, nodeTitle = c("id", "group"))
#' # Map relative abundance to edge color
#' forceGraph(graph_object = mcg, linkStroke = list("relative", "linear", c(0,0.1), c("grey", "blue")))
#' # Set grouping variable
#' forceGraph(graph_object = mcg, groupingVariable = "group")
#' # Disable highlighting of ASV abundance in samples
#' forceGraph(graph_object = mcg, highlightAbundance = FALSE)
#'
#' @import htmlwidgets, crosstalk
#' @export
forceGraph <- function(graph_object,
                       source_color = "type",
                       target_color = "type",
                       groupingVariable = NULL,
                       forceOnGroup = FALSE,
                       nodeStroke = "black",
                       nodeStrokeWidth = 1.5,
                       nodeRadius_source = 10,
                       nodeRadius_target = 4,
                       linkStrokeWidth = 1.5,
                       linkStroke = "gray",
                       linkStrokeOpacity = 0.6,
                       nodeOpacity = 1,
                       nodeStrokeOpacity = 1,
                       nodeTitle = NULL,
                       highlightAbundance = TRUE,
                       width = NULL,
                       height = NULL) {

  # Is crosstalk enabled?
  if (crosstalk::is.SharedData(graph_object$links)) {
    links <- graph_object$links$origData()
    c_key <- graph_object$links$key()
    c_group <- graph_object$links$groupName()
  } else {
    links <- graph_object$links
    c_key <- NULL
    c_group <- NULL
  }

  # Convert to JSON
  g <- list(nodes=graph_object$nodes, links=links)
  g_json <- toJSON(g, force=TRUE)

  # create a list that contains the settings
  settings <- list(
    source_color = source_color,
    target_color = target_color,
    groupingVariable = groupingVariable,
    forceOnGroup = forceOnGroup,
    nodeStroke = nodeStroke,
    nodeRadius_source = nodeRadius_source,
    nodeRadius_target = nodeRadius_target,
    linkStrokeWidth = linkStrokeWidth,
    linkStroke = linkStroke,
    linkStrokeOpacity = linkStrokeOpacity,
    nodeStrokeWidth = nodeStrokeWidth,
    nodeOpacity = nodeOpacity,
    nodeStrokeOpacity = nodeStrokeOpacity,
    nodeTitle = nodeTitle,
    highlightAbundance = highlightAbundance,
    crosstalkKey = c_key,
    crosstalkGroup = c_group
  )

  # pass the data and settings using 'x'
  x <- list(
    graph = g_json,
    settings = settings
  )

  # create the widget
  htmlwidgets::createWidget(
    name = "forceGraph",
    x = x,
    width = width,
    height = height,
    htmlwidgets::sizingPolicy(padding = 10, browser.fill = TRUE),
    dependencies = crosstalk::crosstalkLibs(),
    package = "snowflake"
  )
}

# boilerplate Shiny output and render functions
#' @export
forceGraphOutput <- function(outputId, width = "100%", height = "500px") {
  shinyWidgetOutput(outputId, "forceGraph", width, height, package = "snowflake")
}

#' @export
renderForceGraph <- function(expr, env = parent.frame(), quoted = FALSE) {
  if (!quoted) { expr <- substitute(expr) } # force quoted
  shinyRenderWidget(expr, forceGraphOutput, env, quoted = TRUE)
}
