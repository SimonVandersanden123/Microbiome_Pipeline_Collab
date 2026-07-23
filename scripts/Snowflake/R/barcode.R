#' @title barcode
#' 
#' @description barcode visualization of the microbiome composition graph that shows all observed ASVs/OTUs in the data.
#' 
#' @param graph_object microbiome composition graph build with build_mc_graph().
#' @param color color variable to set color on cells.
#' @param width width of the graph. default is full available width.
#' @param haight hieght of the graph. default is full available height.
#' 
#' @return interactive figure
#' 
#' @examples 
#' barcode(graph_object = mcg, color = "black")
#' 
#' @import htmlwidgets, crosstalk
#' @export
barcode <- function(graph_object,
                       color = NULL,
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
    color = color,
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
    name = "barcode",
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
barcodeOutput <- function(outputId, width = "100%", height = "500px") {
  shinyWidgetOutput(outputId, "barcode", width, height, package = "snowflake")
}

#' @export
renderBarcode <- function(expr, env = parent.frame(), quoted = FALSE) {
  if (!quoted) { expr <- substitute(expr) } # force quoted
  shinyRenderWidget(expr, barcodeOutput, env, quoted = TRUE)
}
