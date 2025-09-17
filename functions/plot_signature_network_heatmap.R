#' Plot Extended Signature Network Heatmap
#'
#' Plots a ComplexHeatmap for an extended gene signature network, with annotation tracks for subtype, signature score, and degree.
#'
#' @param nodes_object List. Output from expand_and_plot_signature_network (must contain node_metrics, signature_genes, expanded_genes, and igraph object 'g').
#' @param expr_data Matrix. Gene expression matrix (genes x samples), rownames are gene symbols, colnames are sample IDs.
#' @param subtype_vector Named character vector. Predicted subtypes for each sample (names must match sample IDs).
#' @param lund_colors Named vector. Colors for each subtype.
#' @param outpath Character. Path to save the PDF output.
#' @param plot_width Numeric. Width of the PDF plot (in inches).
#' @param plot_height Numeric. Height of the PDF plot (in inches).
#' @param plot_title Character. Title for the heatmap plot.
#' @return No return value. The heatmap is saved to the specified outpath.
#' 
#' @import ComplexHeatmap Cairo circlize igraph
#' 
#' @export
#' 
plot_signature_network_heatmap <- function(nodes_object,
                                           expr_data,
                                           subtype_vector,
                                           lund_colors,
                                           cluster_object, 
                                           cluster_colors = c("1" = "#77BEF0", "2" = "#FFCB61", "3" = "#EA5B6F"),
                                           global_signature_cols = FALSE,
                                           top_hubs,
                                           only_hubs = FALSE,
                                           outpath,
                                           plot_width = 20,
                                           plot_height = 20,
                                           plot_title = "Extended Signature Network"){
  
  #check if the number of color match the number of clusters
  if(length(unique(cluster_object)) > length(unique(cluster_colors))) {
    stop(
      paste0(
        "Not enough colors provided for annotating clusters: ",
        length(unique(cluster_object)), " clusters detected, but only ",
        length(unique(cluster_colors)), " colors provided."
      )
    )
  }
  
  #genes to plot
  node_metrics = nodes_object$node_metrics
  
  # Identify rows with any NA in node_metrics
  na_rows <- apply(node_metrics, 1, function(x) any(is.na(x)))
  removed_genes <- rownames(node_metrics)[na_rows]
  
  # Print removed gene names
  if (length(removed_genes) > 0) {
    message("Removed genes with NA in node_metrics: ", paste(removed_genes, collapse = ", "))
    # Remove from expanded_genes as well
    nodes_object$expanded_genes <- setdiff(nodes_object$expanded_genes, removed_genes)
  }
  
  # Remove those rows from node_metrics
  node_metrics <- node_metrics[!na_rows, ]
  
  #deal with only hubs parameter
  if(only_hubs){
    top_hubs = head(node_metrics[order(-node_metrics$degree), "name"], top_hubs)
    genes_to_plot = intersect(top_hubs, rownames(expr_data))
  } else {
    genes_to_plot = nodes_object$expanded_genes
  }
  
  #calculate signature scores for the clusters
  gene_clusters = cluster_object[genes_to_plot]
  
  #or each cluster, calculate the mean expression per sample
  cluster_ids = sort(unique(gene_clusters))
  meta_gene_mat = sapply(cluster_ids, function(cl) {
    genes_in_cl = names(gene_clusters)[gene_clusters == cl]
    if (length(genes_in_cl) == 1) {
      #if only one gene, just use its expression
      expr_data[genes_in_cl, ]
    } else {
      colMeans(expr_data[genes_in_cl, , drop = FALSE])
    }
  })
  
  #ensure meta_gene_mat is samples x clusters
  meta_gene_mat = t(meta_gene_mat)
  meta_gene_mat = t(meta_gene_mat)
  
  #convert to df
  meta_gene_df = as.data.frame(meta_gene_mat)
  
  #sort the cluster IDs
  cluster_ids = sort(unique(gene_clusters))
  
  #dynamically create names like "Cluster 1", "Cluster 2", ...
  colnames(meta_gene_df) = paste0("Cluster ", cluster_ids)
  
  #subtype colors
  col_list = list(Subtype = lund_colors,
                  Cluster = cluster_colors)
  
  #get labels
  pred_lab = subtype_vector
  
  #get split factor
  split = factor(pred_lab, levels = c("UroA", "UroB", "UroC", "GU", "BaSq"))

  #only keep genes present in both cluster_object and genes_to_plot
  if(exists("cluster_object")){
    common_genes = intersect(names(cluster_object), genes_to_plot)
    custom_row_order = names(sort(cluster_object[common_genes], decreasing = FALSE))
    
    #create a factor for cluster assignment in the same order as genes_to_plot
    row_split = factor(cluster_object[genes_to_plot])
    
  }else{
    custom_row_order = genes_to_plot
  }
  
  #heatmap colors
  col_fun = circlize::colorRamp2(c(-2, 0, 2), c("#4DF76F", "black", "#F74D4D"))
  
  #subset expression matrix to genes in network
  expr_subset = expr_data[rownames(expr_data) %in% genes_to_plot, , drop = FALSE]
  expr_data = as.matrix(expr_data)
  
  #calculate signature score (median across selected genes)
  signature_score = apply(expr_subset, 2, median, na.rm = TRUE)
  sig_min = min(signature_score)
  sig_max = max(signature_score)
  col_signature = circlize::colorRamp2(c(sig_min, sig_max), c("#1E93AB", "#E62727"))
  
  #degree annotation
  degree_vec = igraph::degree(nodes_object$g)
  degree_vec = degree_vec[genes_to_plot]
  deg_min = min(degree_vec, na.rm = TRUE)
  deg_max = max(degree_vec, na.rm = TRUE)
  col_degree = circlize::colorRamp2(c(deg_min, deg_max), c("#BCA88D", "#3E3F29"))
  
  if (global_signature_cols) {
    # Global scaling across all clusters
    all_vals <- unlist(meta_gene_df)
    global_min <- min(all_vals, na.rm = TRUE)
    global_max <- max(all_vals, na.rm = TRUE)
    global_mid <- (global_min + global_max) / 2
    
    meta_gene_col_fun <- colorRamp2(
      c(global_min, global_mid, global_max),
      c("#2166AC", "white", "#B2182B")
    )
    
    meta_gene_colors <- setNames(
      rep(list(meta_gene_col_fun), ncol(meta_gene_df)),
      colnames(meta_gene_df)
    )
  } else {
    # Per-cluster scaling
    meta_gene_colors <- setNames(
      lapply(seq_along(meta_gene_df), function(i) {
        vals <- meta_gene_df[[i]]
        rng <- range(vals, na.rm = TRUE)
        mid <- mean(rng)
        colorRamp2(
          c(rng[1], mid, rng[2]),
          c("#2166AC", "white", "#B2182B")
        )
      }),
      colnames(meta_gene_df)
    )
  }
  
  #bind color objects
  col_list = c(
    col_list,
    `Signature Score` = col_signature,
    Degree = col_degree)
  
  #get sample order
  sample_order = get_sample_order(expr_data = expr_data)
  
  #top annotation
  hm_annotation_top = HeatmapAnnotation(
    `Subtype` = pred_lab,
    annotation_name_side = "left",
    gap = unit(1, "mm"),
    simple_anno_size = unit(4, "mm"),
    simple_anno_size_adjust = TRUE,
    col = col_list,
    show_legend = FALSE,
    border = TRUE,
    annotation_name_gp = gpar(fontsize = 8)
  )
  
  #bottom annotation
  bottom_anno = HeatmapAnnotation(
    df = meta_gene_df,
    col = meta_gene_colors,
    annotation_name_side = "left",
    show_legend = FALSE,
    show_annotation_name = TRUE)
  
  #left annotation
  hm_annotation_left = rowAnnotation(
    Degree = degree_vec,
    Cluster = row_split,
    Cluster = row_split,
    col = col_list,
    show_annotation_name = FALSE
  )

  #draw signature heatmap
  hm_signature = Heatmap(
    expr_data[genes_to_plot, , drop = FALSE],
    name = "Expression",
    top_annotation = hm_annotation_top,
    bottom_annotation = bottom_anno,
    left_annotation = hm_annotation_left,
    col = col_fun,
    column_split = split,
    column_title = plot_title,
    column_title_gp = gpar(fontsize = 16, fontface = "bold"),
    cluster_row_slices = FALSE,
    cluster_column_slices = FALSE,
    cluster_columns = FALSE,
    column_order = sample_order,
    height = unit(100, "mm"),
    #width = unit(800, "mm"),
    row_names_side = "left",
    row_title = NULL,
    show_column_names = FALSE,
    show_row_names = TRUE,
    show_row_dend = FALSE,
    border = TRUE,
    row_names_gp = gpar(
      fontsize = 8,
      fontface = ifelse(rownames(expr_data[genes_to_plot, , drop = FALSE]) %in% nodes_object$signature_genes, "bold", "plain"),
      col = ifelse(rownames(expr_data[genes_to_plot, , drop = FALSE]) %in% nodes_object$signature_genes, "#E62727", "black")
    ),
    row_title_gp = gpar(fontsize = 10),
    cluster_rows = FALSE, 
    row_order = custom_row_order,
    row_split = row_split,
    clustering_method_rows = "ward.D2",
    border_gp = gpar(lwd = 1),
    show_heatmap_legend = TRUE,
    row_title_rot = 0
  )
  
  #export plot
  CairoPDF(outpath, width = plot_width, height = plot_height)
  draw(hm_signature)
  dev.off()
}