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
                                           outpath,
                                           plot_width = 20,
                                           plot_height = 20,
                                           plot_title = "Extended Signature Network"){
  
  # Genes to plot
  genes_to_plot <- nodes_object$expanded_genes
  
  # Subtype colors
  col_list <- list(Subtype = lund_colors)
  
  # Get labels
  pred_lab <- subtype_vector
  
  # Get split factor
  split <- factor(pred_lab, levels = c("UroA", "UroB", "UroC", "GU", "BaSq"))
  
  # Heatmap colors
  col_fun <- circlize::colorRamp2(c(-2, 0, 2), c("#4DF76F", "black", "#F74D4D"))
  
  # Subset expression matrix to genes in network
  expr_subset <- expr_data[rownames(expr_data) %in% genes_to_plot, , drop = FALSE]
  expr_data <- as.matrix(expr_data)
  
  # Calculate signature score (median across selected genes)
  signature_score <- apply(expr_subset, 2, median, na.rm = TRUE)
  sig_min <- min(signature_score)
  sig_max <- max(signature_score)
  col_signature <- circlize::colorRamp2(c(sig_min, sig_max), c("#1E93AB", "#E62727"))
  
  # Degree annotation
  degree_vec <- igraph::degree(nodes_object$g)
  degree_vec <- degree_vec[genes_to_plot]
  deg_min <- min(degree_vec, na.rm = TRUE)
  deg_max <- max(degree_vec, na.rm = TRUE)
  col_degree <- circlize::colorRamp2(c(deg_min, deg_max), c("#BCA88D", "#3E3F29"))
  
  # Bind color objects
  col_list <- c(
    col_list,
    `Signature Score` = col_signature,
    Degree = col_degree
  )
  
  # Get sample order
  sample_order <- get_sample_order(expr_data = expr_data)
  
  # Top annotation
  hm_annotation_top <- HeatmapAnnotation(
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
  
  # Bottom annotation
  hm_annotation_bottom <- HeatmapAnnotation(
    `Signature Score` = signature_score,
    annotation_name_side = "left",
    gap = unit(1, "mm"),
    simple_anno_size = unit(4, "mm"),
    simple_anno_size_adjust = TRUE,
    col = col_list,
    show_legend = TRUE,
    border = TRUE,
    annotation_name_gp = gpar(fontsize = 8)
  )
  
  # Left annotation
  hm_annotation_left <- rowAnnotation(
    Degree = degree_vec,
    col = col_list,
    show_annotation_name = FALSE
  )
  
  # Draw signature heatmap
  hm_signature <- Heatmap(
    expr_data[genes_to_plot, , drop = FALSE],
    name = "Expression",
    top_annotation = hm_annotation_top,
    bottom_annotation = hm_annotation_bottom,
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
    cluster_rows = TRUE,
    clustering_distance_rows = "spearman",
    clustering_method_rows = "ward.D2",
    border_gp = gpar(lwd = 1),
    show_heatmap_legend = TRUE,
    row_title_rot = 0
  )
  
  # Export plot
  CairoPDF(outpath, width = plot_width, height = plot_height)
  draw(hm_signature)
  dev.off()
}