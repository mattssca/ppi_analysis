# Create domain-protein network
create_domain_protein_network <- function(subtype_domains, subtype_name) {
  
  library(igraph)
  
  lund_colors = LundTax2023Classifier::lund_colors$lund_colors
  
  # Get domain-protein edges
  edges <- subtype_domains[[subtype_name]] %>%
    dplyr::select(hgnc_symbol, interpro_short_description) %>%
    filter(!is.na(interpro_short_description))
  
  # Create bipartite graph
  g <- graph_from_data_frame(edges, directed = FALSE)
  
  # Set node types
  V(g)$type <- V(g)$name %in% edges$interpro_short_description
  V(g)$node_type <- ifelse(V(g)$type, "domain", "protein")
  
  # Set colors
  V(g)$color <- ifelse(V(g)$type, "lightblue", lund_colors[subtype_name])
  
  return(g)
}

# Create domain networks for each subtype
domain_networks <- list()

for(subtype in names(subtype_domains)) {
  domain_networks[[subtype]] <- create_domain_protein_network(subtype_domains, subtype)
}

# Plot domain network for one subtype
plot_domain_network <- function(domain_graph, subtype_name) {
  
  library(ggraph)
  
  ggraph(domain_graph, layout = "fr") +
    geom_edge_link(alpha = 0.3) +
    geom_node_point(aes(color = node_type), size = 3) +
    geom_node_text(aes(label = name), size = 2, repel = TRUE) +
    scale_color_manual(values = c("domain" = "lightblue", 
                                  "protein" = lund_colors[subtype_name])) +
    labs(title = paste("Domain-Protein Network -", toupper(subtype_name))) +
    theme_graph()
}

# Example: plot domain network for uro subtype
lund_colors = LundTax2023Classifier::lund_colors$lund_colors

plot_domain_network(domain_networks$uro, "uro")
