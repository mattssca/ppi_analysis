# Step 5: Create the visualization function
create_full_ppi_plot <- function(top_n = 100, layout_type = "fr", show_labels = TRUE, highlight_top_hubs = 10) {
  
  # Get top N nodes by degree for cleaner visualization
  degrees <- degree(full_ppi_network)
  top_nodes <- names(sort(degrees, decreasing = TRUE)[1:top_n])
  subgraph <- induced_subgraph(full_ppi_network, top_nodes)
  
  # Identify top hub genes in the subgraph
  subgraph_degrees <- degree(subgraph)
  top_hub_genes <- names(sort(subgraph_degrees, decreasing = TRUE)[1:highlight_top_hubs])
  
  # Convert to tidygraph
  tidy_subgraph <- as_tbl_graph(subgraph) %>%
    activate(nodes) %>%
    mutate(
      degree = degree(subgraph),
      # Get the attributes we added
      dominant_subtype = V(subgraph)$dominant_subtype,
      max_expression = V(subgraph)$max_expression,
      node_color = V(subgraph)$node_color,
      # Add hub highlighting
      is_top_hub = name %in% top_hub_genes
    )
  
  # Create the plot
  p <- ggraph(tidy_subgraph, layout = layout_type) +
    geom_edge_link(alpha = 0.3, color = "gray70", width = 0.5) +
    geom_node_point(aes(size = degree, alpha = max_expression, 
                        stroke = ifelse(is_top_hub, 2, 0.5)), 
                    color = "black", fill = V(subgraph)$node_color, 
                    shape = 21) +
    scale_alpha_continuous(range = c(0.4, 1.0), name = "Max Expression") +
    scale_size_continuous(range = c(3, 12), name = "Degree") +
    labs(title = "Cell Cycle PPI Network - Colored by Dominant Subtype",
         subtitle = paste("Top", top_n, "nodes by degree. Node color = subtype with highest expression. Thick outline =", highlight_top_hubs, "top hubs"),
         caption = "Green=Uro, Blue=GU, Red=BaSq, Orange=Mes, Purple=ScNE") +
    theme_graph() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, family = "Arial"),
      plot.subtitle = element_text(size = 12, hjust = 0.5, family = "Arial"),
      plot.caption = element_text(size = 10, hjust = 0.5, family = "Arial"),
      legend.position = "none"
    )
  
  # Add labels if requested
  if (show_labels) {
    p <- p + geom_node_text(aes(label = name), size = 2.5, repel = TRUE, 
                            max.overlaps = 30, color = "black", fontface = "bold",
                            family = "Arial")
  }
  
  return(p)
}