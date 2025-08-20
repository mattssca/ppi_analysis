#' Plot Cell Cycle PPI Network with Expression-Based Coloring
#'
#' @param ppi_data List object created by create_cell_cycle_ppi()
#' @param subtype Character string specifying the subtype ("uro", "gu", "basq", "mes", "scne")
#' @param n_hubs Integer specifying number of top hub genes to highlight with different fill color
#' @param n_genes Integer specifying number of genes to include in the figure (default 100)
#' @param layout_type Character string specifying layout algorithm ("fr", "kk", "circle", etc.)
#' @param show_labels Logical indicating whether to show gene labels
#' @param save_plot Logical indicating whether to save the plot as PDF
#' @param hub_color Character string specifying color for hub nodes (default "gold")
#' @return ggplot object of the network
#' 
plot_cell_cycle_ppi <- function(ppi_data,
                                subtype = "uro",
                                n_hubs = 10,
                                n_genes = 100,
                                layout_type = "fr",
                                show_labels = TRUE,
                                save_plot = FALSE,
                                hub_color = "gold") {
  
  library(ggraph)
  library(ggplot2)
  library(dplyr)
  library(tidygraph)
  
  # Validate inputs
  if (!subtype %in% names(expr_subtypes)) {
    stop("Invalid subtype. Must be one of: ", paste(names(expr_subtypes), collapse = ", "))
  }
  
  if (n_genes <= 0) {
    stop("n_genes must be a positive integer")
  }
  
  cat("Creating plot for", toupper(subtype), "subtype...\n")
  cat("Including top", n_genes, "genes by degree centrality\n")
  
  # Extract data from ppi_data object
  ppi_graph <- ppi_data$graph
  degrees <- ppi_data$degrees
  hub_scores <- ppi_data$hub_scores
  
  # Get top hubs
  top_hubs <- names(sort(hub_scores))[1:min(n_hubs, length(hub_scores))]
  
  # Get expression data for this subtype
  subtype_expr <- expr_subtypes[[subtype]]
  
  # Calculate mean expression for each gene
  gene_mean_expr <- rowMeans(subtype_expr, na.rm = TRUE)
  
  # Match expression to network nodes
  node_expr <- gene_mean_expr[V(ppi_graph)$gene_symbol]
  names(node_expr) <- V(ppi_graph)$gene_symbol
  
  # Handle missing expression values
  node_expr[is.na(node_expr)] <- median(gene_mean_expr, na.rm = TRUE)
  
  # Get top N genes by degree for visualization
  max_genes <- min(n_genes, vcount(ppi_graph))
  top_degree_nodes <- names(sort(degrees, decreasing = TRUE))[1:max_genes]
  subgraph <- induced_subgraph(ppi_graph, top_degree_nodes)
  
  cat("Displaying top", max_genes, "genes in the network visualization\n")
  
  # Update data for subgraph
  sub_degrees <- degree(subgraph)
  sub_expr <- node_expr[V(subgraph)$gene_symbol]
  sub_hubs <- intersect(top_hubs, V(subgraph)$gene_symbol)
  
  # Create tidygraph object
  tidy_graph <- as_tbl_graph(subgraph) %>%
    activate(nodes) %>%
    mutate(
      gene_symbol = V(subgraph)$gene_symbol,
      degree = sub_degrees,
      mean_expr = sub_expr,
      is_hub = gene_symbol %in% sub_hubs,
      node_fill = ifelse(is_hub, hub_color, "white"),
      expr_category = case_when(
        mean_expr <= quantile(sub_expr, 0.33, na.rm = TRUE) ~ "Low",
        mean_expr <= quantile(sub_expr, 0.67, na.rm = TRUE) ~ "Medium", 
        TRUE ~ "High"
      )
    )
  
  # Create the plot with fill color highlighting
  p <- ggraph(tidy_graph, layout = layout_type) +
    geom_edge_link(alpha = 0.3, color = "gray50", width = 0.5) +
    geom_node_point(aes(size = degree, color = mean_expr, fill = node_fill), 
                    alpha = 0.8, shape = 21, stroke = 1) +
    scale_color_gradient2(low = "green", mid = "black", high = "red",
                          midpoint = median(sub_expr, na.rm = TRUE),
                          name = "Expression\nLevel") +
    scale_fill_identity() +  # Use the actual fill colors specified
    scale_size_continuous(range = c(2, 8), name = "Degree") +
    labs(title = paste("Cell Cycle PPI Network -", toupper(subtype), "Subtype"),
         subtitle = paste("Showing", max_genes, "genes |", n_hubs, "hubs highlighted |", 
                          length(sub_hubs), "hubs visible")) +
    theme_graph() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "gray60"),
      legend.position = "right"
    )
  
  # Add labels if requested
  if (show_labels) {
    p <- p + geom_node_text(aes(label = gene_symbol), size = 2.5, repel = TRUE, 
                            max.overlaps = 20, color = "black", fontface = "bold")
  }
  
  # Save plot if requested
  if (save_plot) {
    filename <- paste0("cell_cycle_ppi_", subtype, "_", max_genes, "genes_", n_hubs, "hubs.pdf")
    ggsave(filename, plot = p, width = 12, height = 10, dpi = 300)
    cat("Plot saved as:", filename, "\n")
  }
  
  # Print hub information
  cat("\nTop", n_hubs, "hub genes in", toupper(subtype), ":\n")
  hub_info <- data.frame(
    Gene = top_hubs[1:min(n_hubs, length(top_hubs))],
    Hub_Score = round(hub_scores[top_hubs[1:min(n_hubs, length(top_hubs))]], 2),
    Expression = round(node_expr[top_hubs[1:min(n_hubs, length(top_hubs))]], 2),
    In_Plot = top_hubs[1:min(n_hubs, length(top_hubs))] %in% sub_hubs,
    stringsAsFactors = FALSE
  )
  print(hub_info)
  
  return(p)
}
