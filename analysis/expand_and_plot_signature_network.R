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
  verbose = TRUE,
  return_data = FALSE,
  string_db = NULL,
  layout_method = "kk",
  show_labels = TRUE,
  lund_colors = FALSE
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
  if (!requireNamespace("viridis", quietly = TRUE)) stop("Please install viridis package.")
  library(viridis)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  if(verbose){
    print(paste("Signature:", signature_name))
  }

  # Get signature genes
  signature_genes <- signature_list[[signature_name]]
  if (is.null(signature_genes)) stop("Signature not found in signature_list.")
  
  if(verbose){
    message("1/8 Running StringDB...")
  }

  # Use provided STRINGdb object or create a new one
  if (is.null(string_db)) {
    string_db <- STRINGdb$new(version="11.5", species=9606, score_threshold=string_score_threshold, input_directory="")
  }
    # Get signature genes
    signature_genes <- signature_list[[signature_name]]
    if (is.null(signature_genes)) stop("Signature not found in signature_list.")

  # Query STRING for neighbors
    if(verbose){
      message("2/8 Mapping genes...")
    }
  mapped <- string_db$map(data.frame(gene=signature_genes), "gene", removeUnmappedRows=TRUE)
  if (nrow(mapped) == 0) stop("None of the signature genes could be mapped to STRING IDs.")
  signature_ids <- mapped$STRING_id
  if (length(signature_ids) == 0 || any(is.na(signature_ids))) stop("No valid STRING IDs for signature genes.")

  if(verbose){
    message("3/8 Finding Neighbours...")
  }
  neighbors <- tryCatch(string_db$get_neighbors(signature_ids), error=function(e) character(0))
  if (length(neighbors) == 0 || all(is.na(neighbors))) stop("No neighbors found for signature genes in STRING.")

  # Only use valid STRING IDs for alias mapping
  valid_neighbors <- neighbors[!is.na(neighbors) & neighbors != ""]
  if (length(valid_neighbors) == 0) stop("No valid STRING neighbor IDs for alias mapping.")

  # Ensure neighbors are character and not NA/empty
  neighbors <- as.character(neighbors)
  neighbors <- neighbors[!is.na(neighbors) & neighbors != ""]
  if (length(neighbors) == 0) stop("No valid STRING neighbor IDs for alias mapping.")

  
  if(verbose){
    message("4/8 Mapping Protein IDs to Gene Symbols...")
  }
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

  
  if(verbose){
    message("5/8 Filtering Edges...")
  }
  
  # Filter by degree (number of connections to signature)
  edges <- suppressWarnings(string_db$get_interactions(signature_ids))
  if (nrow(edges) == 0) stop("No interactions found for signature genes in STRING.")
  neighbor_counts <- table(c(edges$from, edges$to))
  neighbor_counts <- neighbor_counts[names(neighbor_counts) %in% neighbor_genes]
  top_neighbors <- names(sort(neighbor_counts, decreasing=TRUE))[1:max_added_genes]
  expanded_genes <- unique(c(signature_genes, top_neighbors))
  
  # Get STRING IDs for all mapped neighbors and signature genes
  all_genes <- unique(c(signature_genes, neighbor_genes))
  all_mapped <- suppressMessages(suppressWarnings(string_db$map(data.frame(gene=all_genes), "gene", removeUnmappedRows=TRUE)))
  all_ids <- all_mapped$STRING_id
  if (length(all_ids) == 0 || any(is.na(all_ids))) stop("No valid STRING IDs for expanded gene set.")

  
  if(verbose){
    message("6/8 Expanding Network...")
  }
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
  
  print(paste("Using layout method:", layout_method))

  # Compute fixed layout coordinates for the network
  layout_coords <- switch(
    layout_method,
    fr = igraph::layout_with_fr(g),
    kk = igraph::layout_with_kk(g),
    circle = igraph::layout_in_circle(g),
    lgl = igraph::layout_with_lgl(g),
    dh = igraph::layout_with_dh(g),
    igraph::layout_with_fr(g) # default
  )
  layout_coords <- as.matrix(layout_coords)
  if (ncol(layout_coords) != 2) stop("layout_coords must be a matrix with two columns (x and y coordinates)")
  if (nrow(layout_coords) != length(V(g))) stop("Number of rows in layout_coords must match number of nodes in the graph")
  
  if(verbose){
    message("7/8 Visualization")
  }
  # Compute global min/max expression for all expanded genes and all samples
  all_expr <- expr_data[expanded_genes, , drop=FALSE]
  global_expr_range <- range(all_expr, na.rm=TRUE)
  global_mid_expr <- mean(global_expr_range)

  # Subset expression by subtype
  subtypes <- unique(subtype_vector)
  for (subtype in subtypes) {
    if(verbose){
      print(paste("Plotting ", subtype, " Network..."))
    }
    samples <- names(subtype_vector)[subtype_vector == subtype]
    expr_sub <- expr_data[expanded_genes, samples, drop=FALSE]
    mean_expr <- rowMeans(expr_sub, na.rm=TRUE)
    V(g)$expr <- mean_expr[V(g)$name]
    # Plot
    if (!return_data) {
      pdf_file <- file.path(out_dir, paste0(signature_name, "_", subtype, "_network.pdf"))
      CairoPDF(pdf_file, width=20, height=10)
      signature_nodes <- V(g)$name %in% signature_genes
      V(g)$degree <- degree(g)
      p <- ggraph(g, layout = layout_coords) +
        geom_edge_link(color="#4f4f4f", alpha=0.5, width = 0.1) +
        geom_node_point(aes(color=expr, size=degree, shape=signature_nodes), show.legend=TRUE) +
        (if (lund_colors) {
          scale_color_gradient2(low="green", mid="black", high="red", midpoint=global_mid_expr, limits=global_expr_range, na.value="grey80")
        } else {
          scale_color_viridis(option="C", limits=global_expr_range, na.value="grey80")
        }) +
        scale_size_continuous(range=c(4,12)) +
        scale_shape_manual(values=c(16, 17)) +
        (if (show_labels) geom_node_text(aes(label=name), repel=FALSE, size=6, color="black", fontface="bold") else NULL) +
        labs(
          title = paste0(signature_name, " PPI Network (", subtype, ")"),
          subtitle = paste0("Nodes: ", vcount(g), " | Filters: Minimum Degree: ", min_degree, ", String Score Threshold: ", string_score_threshold, 
                            " | Signature genes highlighted | Node size = degree | Color = expression"),
          color = "Expression",
          size = "Degree"
        ) +
        theme_bw() +
        theme(
          axis.title.x = element_blank(), 
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(), 
          axis.ticks = element_blank(),
          plot.title.position = "plot",
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5)
        )
      print(p)
      dev.off()
    }
  }
  if(verbose){
    cat("Signature genes:", length(signature_genes), "\n")
    cat("Neighbors found:", length(neighbors), "\n")
    cat("Mapped neighbor genes:", length(neighbor_genes), "\n")
    cat("Expanded gene set:", length(expanded_genes), "\n")
    cat("Network nodes:", length(V(g)), "\n")
    print(V(g)$name)
    message("8/8 Success!")
  }
  if (return_data) {
    return(list(
      signature_genes = signature_genes,
      neighbors = neighbors,
      neighbor_genes = neighbor_genes,
      expanded_genes = expanded_genes,
      g = g,
      layout_coords = layout_coords
    ))
  } else {
    invisible(NULL)
  }
}
