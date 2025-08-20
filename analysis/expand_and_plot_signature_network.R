#' Expand a gene signature by PPI and plot subtype-specific expression networks
#'
#' @param expr_data Expression matrix (genes x samples)
#' @param subtype_vector Named vector or data frame mapping samples to subtypes
#' @param signature_list List of gene signatures
#' @param signature_name Name of signature to use
#' @param max_added_genes Max number of genes to add to signature (default 20)
#' @param min_degree Minimum number of connections to signature (default 1)
#' @param string_score_threshold Minimum STRINGdb score for added genes (default 400)
#' @param out_dir Output directory for PDFs
#'
#' @return NULL (saves PDFs)
#'
expand_and_plot_signature_network <- function(
  expr_data,
  subtype_vector,
  signature_list,
  signature_name,
  max_added_genes = 20,
  min_degree = 1,
  string_score_threshold = 400,
  out_dir = "viz/signature_networks", 
  verbose = TRUE
) {
  # Load STRINGdb
  if (!requireNamespace("STRINGdb", quietly = TRUE)) stop("Please install STRINGdb package.")
  library(STRINGdb)
  if (!requireNamespace("igraph", quietly = TRUE)) stop("Please install igraph package.")
  library(igraph)
  if (!requireNamespace("ggraph", quietly = TRUE)) stop("Please install ggraph package.")
  library(ggraph)
  if (!requireNamespace("Cairo", quietly = TRUE)) stop("Please install Cairo package.")
  library(Cairo)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # Get signature genes
  signature_genes <- signature_list[[signature_name]]
  if (is.null(signature_genes)) stop("Signature not found in signature_list.")

  # Query STRING for neighbors
  string_db <- STRINGdb$new(version="11.5", species=9606, score_threshold=string_score_threshold, input_directory="")
  mapped <- string_db$map(data.frame(gene=signature_genes), "gene", removeUnmappedRows=TRUE)
  if (nrow(mapped) == 0) stop("None of the signature genes could be mapped to STRING IDs.")
  signature_ids <- mapped$STRING_id
  if (length(signature_ids) == 0 || any(is.na(signature_ids))) stop("No valid STRING IDs for signature genes.")

  neighbors <- tryCatch(string_db$get_neighbors(signature_ids), error=function(e) character(0))
  if (length(neighbors) == 0 || all(is.na(neighbors))) stop("No neighbors found for signature genes in STRING.")

  # Only use valid STRING IDs for alias mapping
  valid_neighbors <- neighbors[!is.na(neighbors) & neighbors != ""]
  if (length(valid_neighbors) == 0) stop("No valid STRING neighbor IDs for alias mapping.")

  # Ensure neighbors are character and not NA/empty
  neighbors <- as.character(neighbors)
  neighbors <- neighbors[!is.na(neighbors) & neighbors != ""]
  if (length(neighbors) == 0) stop("No valid STRING neighbor IDs for alias mapping.")

  # Use biomaRt to map STRING protein IDs to gene symbols
  if (!requireNamespace("biomaRt", quietly = TRUE)) install.packages("biomaRt")
  library(biomaRt)
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  neighbor_peptides <- sub("9606\\.", "", neighbors)
  mapping <- getBM(attributes = c("ensembl_peptide_id", "external_gene_name"),
                   filters = "ensembl_peptide_id",
                   values = neighbor_peptides,
                   mart = mart)
  neighbor_genes <- unique(mapping$external_gene_name)
  neighbor_genes <- setdiff(neighbor_genes, signature_genes)

  # Filter by degree (number of connections to signature)
  edges <- string_db$get_interactions(signature_ids)
  if (nrow(edges) == 0) stop("No interactions found for signature genes in STRING.")
  neighbor_counts <- table(c(edges$from, edges$to))
  neighbor_counts <- neighbor_counts[names(neighbor_counts) %in% neighbor_genes]
  top_neighbors <- names(sort(neighbor_counts, decreasing=TRUE))[1:max_added_genes]
  expanded_genes <- unique(c(signature_genes, top_neighbors))

  # Get STRING IDs for all mapped neighbors and signature genes
  all_genes <- unique(c(signature_genes, neighbor_genes))
  all_mapped <- string_db$map(data.frame(gene=all_genes), "gene", removeUnmappedRows=TRUE)
  all_ids <- all_mapped$STRING_id
  if (length(all_ids) == 0 || any(is.na(all_ids))) stop("No valid STRING IDs for expanded gene set.")

  # Get all interactions among expanded set
  network_edges <- string_db$get_interactions(all_ids)
  if (nrow(network_edges) == 0) stop("No network edges found for expanded gene set.")

  # Count degrees for all nodes in the expanded network
  degree_table <- table(c(network_edges$from, network_edges$to))
  # Map STRING IDs to gene symbols using biomaRt
  node_peptides <- sub("9606\\.", "", names(degree_table))
  node_mapping <- getBM(attributes = c("ensembl_peptide_id", "external_gene_name"),
                        filters = "ensembl_peptide_id",
                        values = node_peptides,
                        mart = mart)
  peptide_to_symbol <- setNames(node_mapping$external_gene_name, node_mapping$ensembl_peptide_id)
  degree_symbols <- peptide_to_symbol[node_peptides]
  # Select top neighbors by degree, excluding signature genes
  top_neighbors <- setdiff(degree_symbols[order(degree_table, decreasing=TRUE)], signature_genes)
  top_neighbors <- unique(top_neighbors[!is.na(top_neighbors)])
  top_neighbors <- head(top_neighbors, max_added_genes)
  expanded_genes <- unique(c(signature_genes, top_neighbors))

  # Build network for expanded genes
  expanded_mapped <- string_db$map(data.frame(gene=expanded_genes), "gene", removeUnmappedRows=TRUE)
  expanded_ids <- expanded_mapped$STRING_id
  network_edges <- string_db$get_interactions(expanded_ids)
  if (nrow(network_edges) == 0) stop("No network edges found for expanded gene set.")
  g <- graph_from_data_frame(network_edges, directed=FALSE)

  # Map node names to gene symbols
  node_peptides <- sub("9606\\.", "", V(g)$name)
  node_mapping <- getBM(attributes = c("ensembl_peptide_id", "external_gene_name"),
                        filters = "ensembl_peptide_id",
                        values = node_peptides,
                        mart = mart)
  peptide_to_symbol <- setNames(node_mapping$external_gene_name, node_mapping$ensembl_peptide_id)
  V(g)$name <- peptide_to_symbol[node_peptides]

  # Compute fixed layout coordinates for the network
  layout_coords <- igraph::layout_with_fr(g)
  # Ensure layout_coords is a matrix with two columns and matches number of nodes
  layout_coords <- as.matrix(layout_coords)
  if (ncol(layout_coords) != 2) stop("layout_coords must be a matrix with two columns (x and y coordinates)")
  if (nrow(layout_coords) != length(V(g))) stop("Number of rows in layout_coords must match number of nodes in the graph")

  print(dim(layout_coords))  # Should be (number of nodes, 2)
  print(layout_coords[1:5, ])  # Should show 5 rows, 2 columns
  print(length(V(g)))         # Should match nrow(layout_coords)
  
  # Subset expression by subtype
  subtypes <- unique(subtype_vector)
  for (subtype in subtypes) {
    samples <- names(subtype_vector)[subtype_vector == subtype]
    expr_sub <- expr_data[expanded_genes, samples, drop=FALSE]
    mean_expr <- rowMeans(expr_sub, na.rm=TRUE)
    V(g)$expr <- mean_expr[V(g)$name]
    # Plot
    pdf_file <- file.path(out_dir, paste0(signature_name, "_", subtype, "_network.pdf"))
    CairoPDF(pdf_file, width=20, height=10)
    # Custom color scale: green (low), black (mid), red (high)
    library(scales)
    expr_range <- range(V(g)$expr, na.rm=TRUE)
    mid_expr <- mean(expr_range)
    p <- ggraph(g, layout = layout_coords) +
      geom_edge_link(color="grey80") +
      geom_node_point(aes(color=expr, size=expr), show.legend=FALSE) +
      scale_color_gradient2(low="green", mid="black", high="red", midpoint=mid_expr, na.value="grey80") +
      scale_size_continuous(range=c(4,12)) +
      geom_node_text(aes(label=name), repel=TRUE, size=3) +
      labs(title=paste0(signature_name, " network: ", subtype)) +
      theme_void()
    print(p)
    dev.off()
  }
  if(verbose){
    cat("Signature genes:", length(signature_genes), "\n")
    cat("Neighbors found:", length(neighbors), "\n")
    cat("Mapped neighbor genes:", length(neighbor_genes), "\n")
    cat("Expanded gene set:", length(expanded_genes), "\n")
    cat("Network nodes:", length(V(g)), "\n")
    print(V(g)$name)
  }
  invisible(NULL)
}
