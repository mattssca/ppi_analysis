#' @title Identify Subtype-Specific Hub Genes Using Network Centrality Analysis
#'
#' @description
#' This function analyzes subtype-specific protein-protein interaction networks 
#' to identify hub genes that are topologically important within each molecular 
#' subtype. The analysis calculates multiple centrality measures and combines 
#' them into a composite hub score to rank genes by their network importance.
#'
#' @details
#' The function uses separate PPI networks for each subtype, where each network 
#' is built from genes characteristic of that subtype. This creates truly different 
#' network topologies per subtype, allowing identification of genes that are both 
#' topologically important and biologically relevant to specific molecular contexts.
#' 
#' The hub score is calculated as the average of three centrality measure ranks:
#' \itemize{
#'   \item Degree centrality: Number of direct protein interactions
#'   \item Betweenness centrality: Frequency of lying on shortest paths between proteins
#'   \item Eigenvector centrality: Connections to other highly connected proteins
#' }
#' 
#' Lower hub scores indicate more central/important genes in the network topology.
#'
#' @param ppi_networks A named list of igraph objects representing protein-protein 
#' interaction networks for each molecular subtype. Each network should have 
#' gene symbols as node names and mean expression values as node attributes 
#' (V(graph)$mean_expr).
#' @param top_n An integer specifying the number of top hub genes to return for 
#' each subtype. Default is 20.
#'
#' @return A named list where each element corresponds to a molecular subtype and 
#'   contains a data.frame with the following columns:
#'   \itemize{
#'     \item gene: Gene symbol
#'     \item degree: Degree centrality value
#'     \item betweenness: Normalized betweenness centrality
#'     \item closeness: Normalized closeness centrality  
#'     \item eigenvector: Eigenvector centrality value
#'     \item pagerank: PageRank centrality value
#'     \item mean_expr: Mean expression value for the gene
#'     \item subtype: Subtype identifier
#'     \item degree_rank: Rank based on degree centrality
#'     \item betweenness_rank: Rank based on betweenness centrality
#'     \item eigenvector_rank: Rank based on eigenvector centrality
#'     \item hub_score: Composite hub score (average of the three centrality ranks)
#'   }
#'   
#'   Results are sorted by hub_score in ascending order (best hubs first).
#'
analyze_subtype_hubs <- function(ppi_networks, 
                                 top_n = 20) {
  # Load required packages
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("igraph package is required but not installed")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr package is required but not installed")
  }
  
  hub_results <- list()
  
  for(subtype in names(ppi_networks)) {
    graph <- ppi_networks[[subtype]]
    
    #calculate multiple centrality measures for this subtype's network
    degrees <- igraph::degree(graph)
    betweenness <- igraph::betweenness(graph, normalized = TRUE)
    closeness <- igraph::closeness(graph, normalized = TRUE)
    eigenvector <- igraph::eigen_centrality(graph)$vector
    pagerank <- igraph::page_rank(graph)$vector
    
    #get mean expression from node attributes
    mean_expr <- igraph::V(graph)$mean_expr
    
    #combine into dataframe - node names are gene symbols
    centrality_df <- data.frame(
      gene = igraph::V(graph)$name,
      degree = degrees,
      betweenness = betweenness,
      closeness = closeness,
      eigenvector = eigenvector,
      pagerank = pagerank,
      mean_expr = mean_expr,
      subtype = subtype
    )
    
    #remove rows with missing gene symbols (if any)
    centrality_df <- centrality_df[!is.na(centrality_df$gene), ]
    
    #rank genes by different centrality measures
    centrality_df$degree_rank <- rank(-centrality_df$degree)
    centrality_df$betweenness_rank <- rank(-centrality_df$betweenness)
    centrality_df$eigenvector_rank <- rank(-centrality_df$eigenvector)
    
    #calculate composite hub score (topology-based)
    centrality_df$hub_score <- (centrality_df$degree_rank + 
                                  centrality_df$betweenness_rank + 
                                  centrality_df$eigenvector_rank) / 3
    
    #get top N hubs
    centrality_df <- centrality_df[order(centrality_df$hub_score), ]
    n_rows <- min(top_n, nrow(centrality_df))
    if (n_rows > 0) {
      centrality_df <- centrality_df[seq_len(n_rows), ]
    }
    
    hub_results[[subtype]] <- centrality_df
    
    #print summary
    cat("Subtype:", subtype, "- Top hub:", centrality_df$gene[1], 
        "with hub score:", round(centrality_df$hub_score[1], 2), "\n")
  }
  
  return(hub_results)
}
