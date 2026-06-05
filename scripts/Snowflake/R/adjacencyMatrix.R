#' @title adjacencyMatrix
#'
#' @description adjacency matrix visualization of the microbiome composition graph
#'
#' @param graph_object microbiome composition graph build with build_mc_graph().
#' @param source_color color variable to set color on row labels (samples).
#' @param target_color color variable to set color on column labels (ASVs).
#' @param cellColor color variable to set color of cells.
#' @param columnBrush add barcode visualization as 2D brush to zoom in on ASVs/OTUs. Boolean.
#' @param width width of the graph. default is full available width.
#' @param haight height of the graph. default is full available height.
#'
#' @return interactive figure
#'
#' @examples
#' # Default
#' adjacencyMatrix(graph_object = mcg)
#' # Color cells by relative abundance
#' adjacencyMatrix(graph_object = mcg, cellColor = "relative")
#' # Color sample labels by group
#' adjacencyMatrix(graph_object = mcg, cellColor = "relative", source_color = "group")
#' # Disable horizontal brushing
#' adjacencyMatrix(graph_object = mcg, cellColor = "relative", columnBrush = FALSE)
#'
#' @import htmlwidgets, crosstalk
#' @export
adjacencyMatrix <- function(graph_object,
                       source_color = NULL,
                       target_color = NULL,
                       cellColor = NULL,
                       columnBrush = TRUE,
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
    cellColor = cellColor,
    columnBrush = columnBrush,
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
    name = "adjacencyMatrix",
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
adjacencyMatrixOutput <- function(outputId, width = "100%", height = "100%") {
  shinyWidgetOutput(outputId, "adjacencyMatrix", width, height, package = "snowflake")
}

#' @export
renderAdjacencyMatrix <- function(expr, env = parent.frame(), quoted = FALSE) {
  if (!quoted) { expr <- substitute(expr) } # force quoted
  shinyRenderWidget(expr, adjacencyMatrixOutput, env, quoted = TRUE)
}
