#' Create Cell Cycle PPI Network Object
#'
#' @param confidence_score Numeric threshold for STRING confidence (default 400)
#' @return List containing igraph object and associated data
#' 
create_cell_cycle_ppi <- function(confidence_score = 400) {
  
  library(STRINGdb)
  library(igraph)
  
  cat("Creating cell cycle PPI network...\n")
  
  # Get STRING database
  string_db <- STRINGdb$new(version = "12.0", species = 9606, score_threshold = confidence_score)
  
  # Map cell cycle genes to STRING IDs
  gene_df <- data.frame(gene = cell_cycle_go, stringsAsFactors = FALSE)
  mapped_genes <- string_db$map(gene_df, "gene", removeUnmappedRows = TRUE)
  
  cat("Mapped", nrow(mapped_genes), "out of", length(cell_cycle_go), "cell cycle genes\n")
  
  # Get PPI network
  ppi_interactions <- string_db$get_interactions(mapped_genes$STRING_id)
  
  # Create igraph network
  ppi_graph <- graph_from_data_frame(ppi_interactions[, c("from", "to")], directed = FALSE)
  
  # Map STRING IDs back to gene symbols
  string_to_gene <- setNames(mapped_genes$gene, mapped_genes$STRING_id)
  V(ppi_graph)$gene_symbol <- string_to_gene[V(ppi_graph)$name]
  
  # Remove nodes without gene symbols
  valid_nodes <- !is.na(V(ppi_graph)$gene_symbol)
  ppi_graph <- induced_subgraph(ppi_graph, which(valid_nodes))
  
  # Calculate centrality measures for the full network
  degrees <- degree(ppi_graph)
  betweenness <- betweenness(ppi_graph)
  eigenvector <- eigen_centrality(ppi_graph)$vector
  
  # Calculate hub scores (lower = better hub)
  degree_rank <- rank(-degrees)
  betweenness_rank <- rank(-betweenness)
  eigenvector_rank <- rank(-eigenvector)
  hub_score <- (degree_rank + betweenness_rank + eigenvector_rank) / 3
  
  cat("PPI network created with", vcount(ppi_graph), "nodes and", ecount(ppi_graph), "edges\n")
  
  # Return list with all necessary data
  ppi_data <- list(
    graph = ppi_graph,
    degrees = degrees,
    betweenness = betweenness,
    eigenvector = eigenvector,
    hub_scores = hub_score,
    gene_symbols = V(ppi_graph)$gene_symbol
  )
  
  return(ppi_data)
}

ppi_network <- create_cell_cycle_ppi(confidence_score = 400)
