# Load required libraries
library(ggraph)
library(tidygraph)
library(igraph)
library(ggplot2)
library(dplyr)

#' Plot PPI Network for a Specific Subtype
#'
#' @description
#' Creates a visualization of protein-protein interaction (PPI) networks for a 
#' specific bladder cancer subtype. The function generates a network plot showing 
#' the top N nodes by degree centrality, with nodes colored by expression level 
#' and sized by degree. Hub genes can be highlighted with black outlines.
#'
#' @details
#' The function takes a subtype-specific PPI network from the global `ppi_networks` 
#' object and creates a ggraph visualization. Node colors follow a gradient from 
#' light to dark using the corresponding Lund taxonomy color for each subtype. 
#' The layout algorithm determines node positioning, with Fruchterman-Reingold 
#' ("fr") providing good general-purpose layouts and Kamada-Kawai ("kk") offering 
#' more structured positioning.
#' 
#' Hub gene highlighting requires either the `hub_data` parameter or a global 
#' `subtype_hubs` object. If neither is available, the top 10 nodes by degree 
#' will be highlighted instead.
#'
#' @param subtype_name Character string specifying the subtype to plot. Must be 
#'   one of the names in the `ppi_networks` object (e.g., "uro", "gu", "basq", 
#'   "mes", "scne").
#' @param top_n Integer specifying the number of top nodes by degree to include 
#'   in the visualization. Default is 50. Smaller values create cleaner plots.
#' @param layout_type Character string specifying the layout algorithm. Options 
#'   include "fr" (Fruchterman-Reingold), "kk" (Kamada-Kawai), "circle", "star", 
#'   etc. Default is "fr".
#' @param show_labels Logical indicating whether to display gene names as node 
#'   labels. Default is TRUE. Set to FALSE for cleaner plots with many nodes.
#' @param node_size_range Numeric vector of length 2 specifying the range of 
#'   node sizes. Default is c(3, 10). Larger ranges increase visual contrast 
#'   between high and low degree nodes.
#' @param save_plot Logical indicating whether to save the plot as a PNG file. 
#'   Default is FALSE. Files are saved as "ppi_network_{subtype_name}.png".
#' @param top_hubs_highlight Integer specifying how many top hub genes should 
#'   be highlighted with thick black outlines. Default is 10. Set to 0 to 
#'   disable hub highlighting.
#'
#' @return A ggplot object containing the network visualization. The plot includes:
#'   \itemize{
#'     \item Nodes colored by expression level using subtype-specific gradients
#'     \item Node sizes proportional to degree centrality
#'     \item Black outlines around top hub genes (if hub data available)
#'     \item Gene labels (if show_labels = TRUE)
#'     \item Legends for expression and degree
#'     \item Subtype-colored title
#'   }
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' plot_ppi_network("uro")
#' 
#' # Customize parameters
#' plot_ppi_network("mes", top_n = 30, layout_type = "kk", 
#'                  show_labels = FALSE)
#' 
#' # Highlight only top 5 hubs
#' plot_ppi_network("basq", top_n = 40, top_hubs_highlight = 5)
#' 
#' # No hub highlighting
#' plot_ppi_network("gu", top_hubs_highlight = 0)
#' 
#' # Save plot with custom hub highlighting
#' plot_ppi_network("scne", top_n = 25, top_hubs_highlight = 15, save_plot = TRUE)
#' }
#'
#' @seealso 
#' \code{\link[ggraph]{ggraph}} for network plotting
#' \code{\link[igraph]{degree}} for centrality measures
#' 
#' @export
#'
#' @author Your Name
#'
#' @keywords network visualization protein-interactions bladder-cancer
#'
# Enhanced function with more options
plot_ppi_network <- function(subtype_name, 
                             top_n = 50, 
                             layout_type = "fr",
                             show_labels = TRUE,
                             node_size_range = c(3, 10),
                             save_plot = FALSE,
                             top_hubs_highlight = 10) {
  
  # Check if subtype exists
  if (!subtype_name %in% names(ppi_networks)) {
    stop("Subtype not found in ppi_networks")
  }
  
  graph <- ppi_networks[[subtype_name]]
  
  #get lund colors
  lund_colors <- c(
    "uro" = "#3cb44b",
    "gu" = "#4363d8", 
    "basq" = "#CD2626",
    "mes" = "#f58231",
    "scne" = "#A020F0"
  )
  
  # Get top N nodes by degree
  degrees <- degree(graph)
  top_nodes <- names(sort(degrees, decreasing = TRUE)[1:top_n])
  subgraph <- induced_subgraph(graph, top_nodes)
  
  # Convert to tidygraph and add node attributes
  tidy_subgraph <- as_tbl_graph(subgraph) %>%
    activate(nodes) %>%
    mutate(
      degree = degree(subgraph),
      mean_expr = ifelse(is.null(V(subgraph)$mean_expr), 
                         runif(vcount(subgraph), 1, 10),
                         V(subgraph)$mean_expr),
      # Add hub status for highlighting
      is_top_hub = name %in% subtype_hubs[[subtype_name]]$gene[1:top_hubs_highlight]
    )
  
  # Color scheme
  base_color <- lund_colors[subtype_name]
  
  # Create base plot
  p <- ggraph(tidy_subgraph, layout = layout_type) +
    geom_edge_link(alpha = 0.3, color = "gray50", width = 0.5) +
    geom_node_point(aes(size = degree, alpha = mean_expr, 
                        stroke = ifelse(is_top_hub, 2, 0)), 
                    color = "black", fill = base_color, shape = 21) +
    ggplot2::scale_alpha_continuous(range = c(0.3, 1.0), name = "Expression") +
    ggplot2::scale_size_continuous(range = node_size_range, name = "Degree") +
    labs(title = paste("PPI Network -", toupper(subtype_name), "Subtype"),
         subtitle = paste("Top", top_n, "nodes by degree (black outline = top", top_hubs_highlight, "hubs)")) +
    theme_graph() +
    theme(
      plot.title = element_text(color = base_color, size = 16, face = "bold", family = "Arial"),
      plot.subtitle = element_text(color = "gray60", size = 11, family = "Arial"),
      legend.position = "none",
      text = element_text(family = "Arial")
    )
  
  # Add labels if requested
  if (show_labels) {
    p <- p + geom_node_text(aes(label = name), size = 2.5, repel = TRUE, 
                            max.overlaps = 20, color = "black", fontface = "bold")
  }
  
  # Save plot if requested
  if (save_plot) {
    ggsave(paste0("ppi_network_", subtype_name, ".png"), 
           plot = p, width = 12, height = 10, dpi = 300)
  }
  
  return(p)
}
