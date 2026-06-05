#' @title build_mcg
#'
#' @description transforms an OTU or ASV data table stored as a phyloseq object into a microbiome composition graph (multivaraite bipartite graph).
#'
#' @param ps_object phyloseq object containing an 'otu_tab' and 'tax_table'.
#' @param beta_dist distance metric to be used to calculate pairwise distances between samples.
#' @param seriation_method seriation algorithm to calculate optimal order of rows and columns in adjacency matrix. See 'seriation' package for available options.
#' @param seriation_dist distance metric to be used to calculate pairwise distances between samples, and to be used as input for the seriation algorithm.
#'
#' @return A list containing two dataframes (nodes and links)
#'
#' @examples
#' mcg <- build_mc_graph(ps_object = ps, beta_dist = "bray", seriation_method = NULL)
#'
#' @export
build_mc_graph <- function(ps_object, beta_dist = "bray", seriation_method = NULL, seriation_dist = "bray") {

  if(ps_object@otu_table@taxa_are_rows) {
    ps_object@otu_table <- t(ps_object@otu_table)
  }

  relAbundance <- transform_sample_counts(ps_object, function(x) x / sum(x) * 100)
  taxonomy <- as.data.frame(ps_object@tax_table)

  beta_dm <- as.matrix(distance(relAbundance, method = beta_dist))
  beta_mds <- cmdscale(beta_dm, k=3)

  columnOrder <- NULL
  rowOrder <- NULL
  if (!is.null(seriation_method)) {
    ASVdistances <- vegdist(t(relAbundance@otu_table), method = seriation_dist)
    columnOrder <- seriate(ASVdistances, method = seriation_method)
    rowOrder <- seriate(as.dist(beta_dm), method = seriation_method)
  }

  # Links
  links <- data.frame("source" = character(), "target" = character(), "absolute" = integer(), "relative" = integer(), stringsAsFactors = FALSE)
  columns <- colnames(ps_object@otu_table)
  rows <- row.names(ps_object@otu_table)
  i <- 1
  for (col in columns) {
    for (row in rows) {
      if (ps_object@otu_table[c(row), c(col)] > 0) {
        links[i, c("source")] <- row
        links[i, c("target")] <- col
        links[i, c("absolute")] <- ps_object@otu_table[c(row), c(col)]
        links[i, c("relative")] <- relAbundance@otu_table[c(row), c(col)]
        links[i, c("index")] <- i
        i <- i + 1
      }
    }
  }
  # Nodes
  node_ids <- unique(c(columns, rows))
  nodes <- data.frame("id" = character(),
                      "type" = character(),
                      "domain" = character(),
                      "phylum" = character(),
                      "class" = character(),
                      "order" = character(),
                      "family" = character(),
                      "genus" = character(),
                      "species" = character(),
                      "index" = character(),
                      "abundance" = integer(),
                      "alphaObserved" = integer(),
                      "alphaChao1" = integer(),
                      "alphaShannon" = integer(),
                      "alphaInvSimpson" = integer(),
                      "observedIn" = integer(),
                      "beta_mds_1" = integer(),
                      "beta_mds_2" = integer(),
                      "beta_mds_3" = integer(),
                      "seriation" = integer(),
                      "connectedWith" = character(),
                      stringsAsFactors = FALSE)
  i <- 1
  for (node in rows) {
    nodes[i, "id"] <- node
    nodes[i, "type"] <- "source"
    nodes[i, "alphaObserved"] <- as.numeric(estimate_richness(ps_object@otu_table[node,], split = TRUE, measures = c("Observed")))
    nodes[i, "alphaChao1"] <- as.numeric(estimate_richness(ps_object@otu_table[node,], measures=c("Chao1"))[1])
    nodes[i, "alphaShannon"] <- as.numeric(estimate_richness(ps_object@otu_table[node,], measures=c("Shannon")))
    nodes[i, "alphaInvSimpson"] <- as.numeric(estimate_richness(ps_object@otu_table[node,], measures=c("InvSimpson")))
    nodes[i, "beta_mds_1"] <- as.numeric(beta_mds[node,1])
    nodes[i, "beta_mds_2"] <- as.numeric(beta_mds[node,2])
    nodes[i, "beta_mds_3"] <- as.numeric(beta_mds[node,3])
    nodes[i, "seriation"] <- if(!is.null(rowOrder)) rowOrder[[1]]$order[i] else NA
    nodes[i, "index"] <- i
    nodes[i, "connectedWith"] <- paste0(links[which(links$source == node), "target"], collapse = ",")
    if(!is.null(ps_object@sam_data)) {
      for(col in names(ps_object@sam_data)) {
        nodes[i, col] <- if(node %in% rownames(ps_object@sam_data)) ps_object@sam_data[which(rownames(ps_object@sam_data) == node), col] else NA
      }
    }
    i <- i + 1
  }
  for (node in columns) {
    j <- 1
    nodes[i, "id"] <- node
    nodes[i, "type"] <- "target"
    nodes[i, "domain"] <- taxonomy[node, "domain"]
    nodes[i, "phylum"] <- taxonomy[node, "phylum"]
    nodes[i, "class"] <- taxonomy[node, "class"]
    nodes[i, "order"] <- taxonomy[node, "order"]
    nodes[i, "family"] <- taxonomy[node, "family"]
    nodes[i, "genus"] <- taxonomy[node, "genus"]
    nodes[i, "species"] <- taxonomy[node, "species"]
    nodes[i, "abundance"] <- colSums(ps_object@otu_table[,node])/sum(as.matrix(ps_object@otu_table))
    nodes[i, "observedIn"] <- length(links$source[which(links$target == node)])
    nodes[i, "seriation"] <- if(!is.null(columnOrder)) columnOrder[[1]]$order[j] else NA
    nodes[i, "connectedWith"] <- paste0(links[which(links$target == node), "source"], collapse = ",")
    nodes[i, "index"] <- i
    j <- j +1
    i <- i + 1
  }

  graph <- list(nodes = nodes, links = crosstalk::SharedData$new(links))
  # graph <- list(nodes = nodes, links = links)
  return(graph)
}
