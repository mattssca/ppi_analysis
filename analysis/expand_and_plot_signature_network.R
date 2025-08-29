#' Expand a gene signature by PPI and plot subtype-specific expression networks
#'
#' This function expands a gene signature using STRINGdb, builds a PPI network, and visualizes
#' subtype-specific expression profiles for the network nodes. Node color reflects mean/median/zscore
#' expression for each gene in the selected subtype, with options for global or per-subtype color scaling.
#' Node size can be fixed or scaled by degree. Returns network metrics if requested.
#'
#' @param expr_data Expression matrix (genes x samples).
#' @param subtype_vector Named vector mapping samples to subtypes.
#' @param signature_list List of gene signatures.
#' @param signature_name Name of the signature to use.
#' @param max_added_genes Maximum number of genes to add to the signature (default: 20).
#' @param min_degree Minimum number of connections to the signature (default: 1).
#' @param string_score_threshold Minimum STRINGdb score for added genes (default: 400).
#' @param out_dir Output directory for saving PDFs (default: "viz/signature_networks").
#' @param verbose Logical; if TRUE, print progress messages (default: TRUE).
#' @param return_data Logical; if TRUE, return network and node metrics instead of plotting (default: FALSE).
#' @param layout_method Network layout method (e.g., "kk", "fr", "circle").
#' @param show_labels Logical; if TRUE, show gene symbol labels on nodes (default: TRUE).
#' @param expr_summary Method for summarizing expression: "mean", "median", or "zscore" (default: "mean").
#' @param node_color Node color scale: "lund", "magma", "scico", "viridis" (default: "lund").
#' @param theme Plot theme: "light" or "dark" (default: "light").
#' @param node_size Fixed node size (default: 10).
#' @param node_degree Logical; if TRUE, scale node size by degree (default: FALSE).
#' @param plot_width Width of the output plot in inches (default: 15).
#' @param plot_height Height of the output plot in inches (default: 10).
#' @param color_scale_global Logical; if TRUE, use a global color scale for all subtypes; if FALSE, scale per subtype (default: TRUE).
#' @param genes_blacklist A vector of gene symbols to exclude from the top neighbors (default: NULL).
#'
#' @return NULL (plots saved to PDF) or a list with network and node metrics if return_data = TRUE.
#'
#' @import scico STRINGdb igraph ggraph Cairo viridis biomaRt
#' @importFrom STRINGdb STRINGdb
#' @importFrom biomaRt useMart getBM
#' @importFrom igraph graph_from_data_frame V V<-
#'
#' @export
#'
expand_and_plot_signature_network <- function(expr_data,
                                              subtype_vector,
                                              signature_list,
                                              signature_name,
                                              max_added_genes = 20,
                                              min_degree = 1,
                                              string_score_threshold = 400,
                                              out_dir = "viz/signature_networks", 
                                              verbose = TRUE,
                                              return_data = FALSE,
                                              layout_method = "kk",
                                              show_labels = TRUE,
                                              expr_summary = "mean",
                                              node_color = "lund",
                                              theme = "light",
                                              node_size = 10,
                                              node_degree = FALSE,
                                              plot_width = 15,
                                              plot_height = 10,
                                              color_scale_global = TRUE,
                                              genes_blacklist = NULL) {

  #suppress unused variable warnings
  all_ids <- NULL
  
  #create directory, if non-existing
  if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  #verbose message
  if(verbose){print(paste("Signature:", signature_name))}
  
  #get signature genes
  signature_genes <- signature_list[[signature_name]]
  if(is.null(signature_genes)) stop("Signature not found in signature_list.")
  
  
  ########### Step 1 ############
  if(verbose){message("1/8 Running StringDB...")}
  
  #run stringDB for retreiving PPI
  string_db <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = string_score_threshold, input_directory = "")
  
  
  ########### Step 2 ############
  if(verbose){message("2/8 Mapping genes...")}
  
  #get STRING IDs for the genes in the selected signature
  mapped <- string_db$map(data.frame(gene = signature_genes), "gene", removeUnmappedRows = TRUE)
  signature_ids <- mapped$STRING_id
  
  
  ########### Step 3 ############
  if(verbose){message("3/8 Finding Neighbours...")}
  
  #get neighbors for the signature genes
  neighbors <- tryCatch(string_db$get_neighbors(signature_ids), error = function(e) character(0))
  if(length(neighbors) == 0 || all(is.na(neighbors))) stop("No neighbors found for signature genes in STRING.")
  
  #ensure neighbors are character and not NA/empty
  neighbors <- as.character(neighbors)
  neighbors <- neighbors[!is.na(neighbors) & neighbors != ""]
  if(length(neighbors) == 0) stop("No valid STRING neighbor IDs for alias mapping.")
  
  
  ########### Step 4 ############
  if(verbose){message("4/8 Mapping Protein IDs to Gene Symbols...")}
  
  #use biomaRt to map STRING protein IDs to gene symbols
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  neighbor_peptides <- sub("9606\\.", "", neighbors)
  
  mapping <- getBM(attributes = c("ensembl_peptide_id", "external_gene_name"),
                   filters = "ensembl_peptide_id",
                   values = neighbor_peptides,
                   mart = mart)
  
  #get unique neighbor genes that are not in the signature object
  neighbor_genes <- unique(mapping$external_gene_name)
  neighbor_genes <- setdiff(neighbor_genes, signature_genes)
  
  
  ########### Step 5 ############
  if(verbose){message("5/8 Filtering Edges...")}
  
  #build initial network (signature genes)
  #filter by degree (number of connections to signature)
  edges <- suppressWarnings(string_db$get_interactions(signature_ids))
  
  if(nrow(edges) == 0) stop("No interactions found for signature genes in STRING.")
  
  #count the number of connections (degree) for each neighbor gene, 
  neighbor_counts <- table(c(edges$from, edges$to))
  
  #filter to include only valid neighbor genes, select the top neighbors 
  neighbor_counts <- neighbor_counts[names(neighbor_counts) %in% neighbor_genes]
  
  #exclude genes in the blacklist
  if(!is.null(genes_blacklist)){
    neighbor_counts <- neighbor_counts[!names(neighbor_counts) %in% genes_blacklist]
  }
  
  #based on degree (up to max_added_genes), and combine them with the 
  top_neighbors <- names(sort(neighbor_counts, decreasing=TRUE))[1:max_added_genes]
  
  #original signature genes to create the expanded gene set.
  expanded_genes <- unique(c(signature_genes, top_neighbors))
  
  #combine signature genes and neighbor genes, map them to STRING IDs, 
  all_genes <- unique(c(signature_genes, neighbor_genes))
  all_mapped <- suppressMessages(suppressWarnings(string_db$map(data.frame(gene = all_genes), "gene", removeUnmappedRows = TRUE)))
  all_ids <- all_mapped$STRING_id
  if(length(all_ids) == 0 || any(is.na(all_ids))) stop("No valid STRING IDs for expanded gene set.")
  
  
  ########### Step 6 ############
  if(verbose){message("6/8 Expanding Network...")}
  
  #get all interactions among expanded set
  network_edges <- string_db$get_interactions(all_ids)
  
  if(nrow(network_edges) == 0) stop("No network edges found for expanded gene set.")
  
  #count degrees for all nodes in the expanded network
  degree_table <- table(c(network_edges$from, network_edges$to))
  
  #map STRING IDs to gene symbols using biomaRt
  node_peptides <- sub("9606\\.", "", names(degree_table))
  
  node_mapping <- getBM(attributes = c("ensembl_peptide_id", "external_gene_name"),
                        filters = "ensembl_peptide_id",
                        values = node_peptides,
                        mart = mart)
  
  peptide_to_symbol <- setNames(node_mapping$external_gene_name, node_mapping$ensembl_peptide_id)
  
  degree_symbols <- peptide_to_symbol[node_peptides]
  
  #select top neighbors by degree, excluding signature genes
  top_neighbors <- setdiff(degree_symbols[order(degree_table, decreasing = TRUE)], signature_genes)
  top_neighbors <- unique(top_neighbors[!is.na(top_neighbors)])
  top_neighbors <- head(top_neighbors, max_added_genes)
  expanded_genes <- unique(c(signature_genes, top_neighbors))
  
  #build network for expanded genes
  expanded_mapped <- string_db$map(data.frame(gene = expanded_genes), "gene", removeUnmappedRows = TRUE)
  expanded_ids <- expanded_mapped$STRING_id
  network_edges <- string_db$get_interactions(expanded_ids)
  
  if(nrow(network_edges) == 0) stop("No network edges found for expanded gene set.")
  
  g <- graph_from_data_frame(network_edges, directed = FALSE)
  
  #map node names to gene symbols
  node_peptides <- sub("9606\\.", "", V(g)$name)
  
  node_mapping <- getBM(attributes = c("ensembl_peptide_id", "external_gene_name"),
                        filters = "ensembl_peptide_id",
                        values = node_peptides,
                        mart = mart)
  
  peptide_to_symbol <- setNames(node_mapping$external_gene_name, node_mapping$ensembl_peptide_id)
  
  V(g)$name <- peptide_to_symbol[node_peptides]
  
  print(paste("Using layout method:", layout_method))
  
  #compute fixed layout coordinates for the network
  layout_coords <- switch(layout_method,
                          fr = igraph::layout_with_fr(g),
                          kk = igraph::layout_with_kk(g),
                          circle = igraph::layout_in_circle(g),
                          lgl = igraph::layout_with_lgl(g),
                          dh = igraph::layout_with_dh(g),
                          igraph::layout_with_fr(g))
  
  layout_coords <- as.matrix(layout_coords)
  if(ncol(layout_coords) != 2) stop("layout_coords must be a matrix with two columns (x and y coordinates)")
  if(nrow(layout_coords) != length(V(g))) stop("Number of rows in layout_coords must match number of nodes in the graph")
  
  
  ########### Step 7 ############
  if(verbose){message("7/8 Visualization")}
  
  #helper function to calculate expression values based on the summary method
  calculate_expression <- function(expr_data, method) {
    if (method == "mean") {
      return(rowMeans(expr_data, na.rm = TRUE))
    } else if (method == "median") {
      return(apply(expr_data, 1, median, na.rm = TRUE))
    } else if (method == "zscore") {
      return(apply(expr_data, 1, function(x) {
        m <- mean(x, na.rm = TRUE)
        s <- sd(x, na.rm = TRUE)
        if (is.na(s) || s == 0) return(0)
        mean((x - m) / s, na.rm = TRUE)
      }))
    } else {
      stop("Unknown expr_summary option. Use 'mean', 'median', or 'zscore'.")
    }
  }
  
  #calculate global expression range
  all_node_expr <- unlist(lapply(unique(subtype_vector), function(subtype) {
    samples <- names(subtype_vector)[subtype_vector == subtype]
    expr_sub <- expr_data[expanded_genes, samples, drop = FALSE]
    calculate_expression(expr_sub, expr_summary)
  }))
  
  global_expr_range <- range(all_node_expr, na.rm = TRUE)
  global_mid_expr <- mean(global_expr_range)
  
  #subset expression by subtype and plot
  for (subtype in unique(subtype_vector)) {
    if (verbose) {
      print(paste("Plotting", subtype, "Network..."))
    }
    
    samples <- names(subtype_vector)[subtype_vector == subtype]
    valid_samples <- intersect(samples, colnames(expr_data))
    valid_genes <- intersect(expanded_genes, rownames(expr_data))
    
    if (length(valid_samples) == 0 || length(valid_genes) == 0) {
      warning(paste("Skipping subtype", subtype, "due to missing data."))
      next
    }
    
    expr_sub <- expr_data[valid_genes, valid_samples, drop = FALSE]
    expr_vals <- calculate_expression(expr_sub, expr_summary)
    
    #map expression values to network nodes
    expr_vec <- setNames(rep(NA, length(V(g)$name)), V(g)$name)
    common_names <- intersect(names(expr_vec), names(expr_vals))
    expr_vec[common_names] <- expr_vals[common_names]
    V(g)$expr <- expr_vec
    
    if (!return_data) {
      pdf_file <- file.path(out_dir, paste0(signature_name, "_", subtype, "_network.pdf"))
      CairoPDF(pdf_file, width=plot_width, height=plot_height)
      signature_nodes <- V(g)$name %in% signature_genes
      V(g)$degree <- degree(g)
      
      if (node_degree) {
        node_size_vec <- V(g)$degree
        size_label <- "Degree"
        size_scale <- scale_size_continuous(range = c(7, 16))
      } else {
        node_size_val <- node_size
        node_size_vec <- rep(node_size_val, length(V(g)))
        size_label <- "Node Size"
        size_scale <- scale_size_continuous(range = c(node_size, node_size))
      }
      
      #plot
      p <- ggraph(g, layout = layout_coords) +
        geom_edge_link(color = "#4f4f4f", alpha = 0.5, width = 0.1) +
        geom_node_point(aes(color = expr, size = node_size_vec, shape = signature_nodes), show.legend = FALSE) +
        (
          if (node_color == "lund") {
            if (color_scale_global) {
              scale_color_gradientn(colors = c("green", "black", "red"),
                                    values = scales::rescale(c(global_expr_range[1], 0, global_expr_range[2])),
                                    limits = global_expr_range,
                                    na.value = "grey80")
            } else {
              subtype_range <- range(expr_vals, na.rm=TRUE)
              scale_color_gradientn(colors = c("green", "black", "red"),
                                    values = scales::rescale(c(subtype_range[1], 0, subtype_range[2])),
                                    limits = subtype_range,
                                    na.value = "grey80")
            }
          } else if (node_color == "magma") {
            scale_color_viridis(option="B", limits=global_expr_range, na.value="grey80")
          } else if (node_color == "scico") {
            scico::scale_color_scico(palette = "roma", limits=global_expr_range, na.value="grey80")
          } else {
            scale_color_viridis(option="E", limits=global_expr_range, na.value="grey80")
          }
        ) +
        size_scale +
        scale_shape_manual(values=c(16, 17)) +
        (
          if (show_labels && theme == "dark") {
            geom_label(aes(x = x, y = y, label = name), fill = "white", color = "black", size = 6, fontface = "bold", label.size = 0.2)
          } else if (show_labels) {
            geom_label(aes(x = x, y = y, label = name), fill = "white", color = "black", size = 6, fontface = "bold", label.size = 0.2)
          } else {
            NULL
          }
        ) +
        labs(
          title = paste0(signature_name, " PPI Network (", subtype, ")"),
          subtitle = paste0("Nodes: ", vcount(g), " | Filters: Minimum Degree: ", min_degree, ", String Score Threshold: ", string_score_threshold, 
                            " | Signature genes highlighted | Node size = ", size_label, " | Color = expression"),
          color = "Expression",
          size = size_label
        ) +
        (
          if (theme == "dark") {
            theme_bw(base_size = 15) +
              theme(
                plot.background = element_rect(fill = "grey20", color = NA),
                panel.background = element_rect(fill = "grey20", color = NA),
                legend.background = element_rect(fill = "grey20", color = NA),
                legend.text = element_text(color = "white"),
                legend.title = element_text(color = "white"),
                axis.title.x = element_blank(), 
                axis.title.y = element_blank(),
                axis.text.x = element_blank(),
                axis.text.y = element_blank(), 
                axis.ticks = element_blank(),
                plot.title.position = "plot",
                plot.title = element_text(hjust = 0.5, face = "bold", color = "white"),
                plot.subtitle = element_text(hjust = 0.5, color = "white")
              )
          } else {
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
          }
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
    message("Compiling data for return...")
    
    #compute additional node metrics
    n_nodes <- length(V(g))
    node_metrics <- data.frame(
      name = V(g)$name,
      degree = degree(g),
      is_signature = V(g)$name %in% signature_genes,
      betweenness = if(length(betweenness(g)) == n_nodes) betweenness(g) else rep(NA, n_nodes),
      closeness = if(length(closeness(g)) == n_nodes) closeness(g) else rep(NA, n_nodes),
      eigenvector = if(length(eigen_centrality(g)$vector) == n_nodes) eigen_centrality(g)$vector else rep(NA, n_nodes),
      hub_score = if(length(hub_score(g)$vector) == n_nodes) hub_score(g)$vector else rep(NA, n_nodes),
      authority_score = if(length(authority_score(g)$vector) == n_nodes) authority_score(g)$vector else rep(NA, n_nodes),
      community = if(length(membership(cluster_louvain(g))) == n_nodes) as.numeric(membership(cluster_louvain(g))) else rep(NA, n_nodes)
    )

    # Add mean expression for each subtype as new columns
    subtypes <- unique(subtype_vector)
    for (subtype in subtypes) {
      samples <- names(subtype_vector)[subtype_vector == subtype]
      valid_samples <- intersect(samples, colnames(expr_data))
      expr_sub <- expr_data[node_metrics$name, valid_samples, drop = FALSE]
      node_metrics[[paste0("mean_expr_", subtype)]] <- rowMeans(expr_sub, na.rm = TRUE)
    }
    
    return(list(
      signature_genes = signature_genes,
      neighbors = neighbors,
      neighbor_genes = neighbor_genes,
      expanded_genes = expanded_genes,
      g = g,
      layout_coords = layout_coords,
      node_metrics = node_metrics
    ))
  } else {
    invisible(NULL)
  }
}
