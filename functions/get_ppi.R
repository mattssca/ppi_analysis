#' @title Create Protein-Protein Interaction Networks from Expression Data using STRING Database
#'
#' @details This function constructs protein-protein interaction (PPI) networks by mapping 
#' gene symbols to the STRING database, retrieving high-confidence interactions, 
#' and building an igraph network object annotated with gene expression values. 
#' The resulting network uses gene symbols as node identifiers for downstream 
#' analysis and visualization.
#'
#' @description
#' The function implements a multi-step process to create PPI networks:
#' 1. Map input gene symbols to STRING database protein identifiers
#' 2. Retrieve protein-protein interactions with confidence scores ≥ 400
#' 3. Convert STRING IDs back to gene symbols for interpretability
#' 4. Build an undirected igraph network object
#' 5. Annotate nodes with mean expression values from input data
#' 6. Optionally return raw interaction data for external visualization
#' 
#' The function uses STRING database version 11.5 with a confidence threshold of 
#' 400 (medium-high confidence) to balance network coverage with interaction reliability.
#'
#' @param expr_data A numeric matrix or data.frame containing gene expression data. 
#' Genes must be in rows with gene symbols as row names, and samples in columns. 
#' Missing values are handled with na.rm = TRUE in mean calculations.
#' @param return_cytoscape Logical value indicating the return format. If TRUE, 
#' returns a data.frame of interactions suitable for Cytoscape import. If FALSE 
#' (default), returns an igraph object for network analysis in R.
#'
#' @return An igraph object with gene symbols as node names.vNode attribute 'mean_expr' 
#' containing mean expression values. Edge attributes from STRING database (confidence scores, 
#' interaction types). Undirected network structure.
#'
#' @references
#' Szklarczyk D, et al. (2021). The STRING database in 2021: customizable 
#' protein-protein networks, and functional characterization of user-uploaded 
#' gene/measurement sets. Nucleic Acids Research, 49(D1), D605-D612.
#'
run_stringdb <- function(expr_data = NULL, 
                         return_cytoscape = FALSE){
  
  # Load required packages
  if (!requireNamespace("STRINGdb", quietly = TRUE)) {
    stop("STRINGdb package is required but not installed")
  }
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("igraph package is required but not installed")
  }
  
  #initialize STRINGdb (for human, species = 9606)
  string_db <- STRINGdb::STRINGdb$new(version = "11.5", species = 9606, score_threshold = 400)
  
  #map gene symbols to STRING IDs (ignore the misleading warning)
  genes <- rownames(expr_data)
  cat("Mapping", length(genes), "genes to STRING database...\n")
  
  #suppress the misleading warning
  mapped <- suppressWarnings(
    string_db$map(data.frame(gene = genes), "gene", removeUnmappedRows = TRUE)
  )
  
  cat("Successfully mapped", nrow(mapped), "out of", length(genes), "genes\n")
  cat("Mapping success rate:", round(nrow(mapped)/length(genes) * 100, 1), "%\n")
  
  #get interactions for mapped STRING IDs
  interactions <- string_db$get_interactions(mapped$STRING_id)
  
  #create lookup: STRING ID -> gene symbol
  string_to_symbol <- setNames(mapped$gene, mapped$STRING_id)
  
  if(return_cytoscape){
    # Add gene symbols to the interactions data frame
    interactions$from_gene <- string_to_symbol[interactions$from]
    interactions$to_gene   <- string_to_symbol[interactions$to]
    return(interactions)
  }
  
  # Convert STRING IDs to gene symbols in interactions
  interactions$from_gene <- string_to_symbol[interactions$from]
  interactions$to_gene <- string_to_symbol[interactions$to]
  
  # Create igraph object using gene symbols instead of STRING IDs
  graph_obj <- igraph::graph_from_data_frame(
    interactions[, c("from_gene", "to_gene", setdiff(names(interactions), c("from", "to", "from_gene", "to_gene")))], 
    directed = FALSE
  )
  
  #calculate mean expression for each gene
  mean_expr <- rowMeans(expr_data, na.rm = TRUE)
  
  # Add expression values (node names are now gene symbols)
  igraph::V(graph_obj)$mean_expr <- mean_expr[igraph::V(graph_obj)$name]
  
  cat("Created network with", igraph::vcount(graph_obj), "nodes and", igraph::ecount(graph_obj), "edges\n")
  
  return(graph_obj)
}
