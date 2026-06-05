forceGraph_SV <- function(graph_object,
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
