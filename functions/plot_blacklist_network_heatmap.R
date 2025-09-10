plot_blacklist_network_heatmap <- function(nodes_object,
                                           expr_data,
                                           subtype_vector,
                                           lund_colors,
                                           subtype_class = "7_class",
                                           blacklist_genes = NULL,
                                           cut_tree = 2,
                                           outpath,
                                           plot_width = 20,
                                           plot_height = 20,
                                           plot_title = "Extended Signature Network"){
  
  #genes to plot
  node_metrics = nodes_object$node_metrics
  
  #genes to plot
  genes_to_plot <- setdiff(nodes_object$expanded_genes, blacklist_genes)

  #subtype colors
  col_list = list(Subtype = lund_colors)
  
  #get labels
  pred_lab = subtype_vector
  
  #get split factor
  if(subtype_class == "5_class"){
    split = factor(pred_lab, levels = c("Uro", "GU", "BaSq"))
    
  }else if(subtype_class == "7_class"){
    split = factor(pred_lab, levels = c("UroA", "UroB", "UroC", "GU", "BaSq"))
  }
  
  #heatmap colors
  col_fun = circlize::colorRamp2(c(-2, 0, 2), c("#4DF76F", "black", "#F74D4D"))
  
  #subset expression matrix to genes in network
  expr_subset = as.matrix(expr_data[rownames(expr_data) %in% genes_to_plot, , drop = FALSE])
  expr_data = as.matrix(expr_data)

  #degree annotation
  degree_vec = igraph::degree(nodes_object$g)
  degree_vec = degree_vec[rownames(expr_subset)]
  deg_min = min(degree_vec, na.rm = TRUE)
  deg_max = max(degree_vec, na.rm = TRUE)
  col_degree = circlize::colorRamp2(c(deg_min, deg_max), c("#BCA88D", "#3E3F29"))

  #bind color objects
  col_list = c(col_list,
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
  
  #left annotation
  hm_annotation_left = rowAnnotation(
    Degree = degree_vec,
    col = col_list,
    show_annotation_name = FALSE
  )

  #temporary draw heatmap to get cluster info
  pdf(NULL)
  hm_signature_tmp = Heatmap(
    expr_subset,
    cluster_rows = TRUE,
    row_split = cut_tree,
    clustering_distance_rows = "pearson",
    clustering_method_rows = "ward.D2"
  )
  
  hm_drawn <- draw(hm_signature_tmp)
  dev.off()
  
  #extract cluster assignments for rows (genes)
  row_clusters <- ComplexHeatmap::row_order(hm_drawn)
  
  #create a named vector mapping gene names to cluster number
  gene_cluster_assignment <- rep(NA, nrow(expr_subset))
  names(gene_cluster_assignment) <- rownames(expr_subset)
  for (i in seq_along(row_clusters)) {
    gene_cluster_assignment[row_clusters[[i]]] <- i
  }
  
  #for each cluster, calculate mean expression per sample
  cluster_ids <- sort(unique(gene_cluster_assignment))
  meta_gene_mat <- sapply(cluster_ids, function(cl) {
    genes_in_cl <- names(gene_cluster_assignment)[gene_cluster_assignment == cl]
    colMeans(expr_subset[genes_in_cl, , drop = FALSE], na.rm = TRUE)
  })
  meta_gene_df <- as.data.frame(meta_gene_mat)
  colnames(meta_gene_df) <- paste0("Cluster ", cluster_ids)
  
  #create color mapping for the meta-gene tracks
  all_vals <- unlist(meta_gene_df)
  min_val <- min(all_vals, na.rm = TRUE)
  max_val <- max(all_vals, na.rm = TRUE)
  mid_val <- (min_val + max_val) / 2
  
  meta_gene_col_fun <- circlize::colorRamp2(
    c(min_val, mid_val, max_val),
    c("#2166AC", "white", "#B2182B")
  )
  meta_gene_colors <- setNames(
    rep(list(meta_gene_col_fun), ncol(meta_gene_df)),
    colnames(meta_gene_df)
  )
  
  #create bottom annotation
  bottom_anno <- HeatmapAnnotation(
    df = meta_gene_df,
    col = meta_gene_colors,
    annotation_name_side = "left",
    show_annotation_name = TRUE, 
    show_legend = FALSE
  )
  
  #redraw the ehatmap with bottom annotation
  hm_signature = Heatmap(
    expr_subset,
    name = "Expression",
    top_annotation = hm_annotation_top, 
    left_annotation = hm_annotation_left, 
    bottom_annotation = bottom_anno,
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
      fontface = ifelse(rownames(expr_subset) %in% nodes_object$signature_genes, "bold", "plain"),
      col = ifelse(rownames(expr_subset) %in% nodes_object$signature_genes, "#E62727", "black")),
    row_title_gp = gpar(fontsize = 10),
    cluster_rows = TRUE, 
    row_split = cut_tree,
    clustering_distance_rows = "pearson",
    clustering_method_rows = "ward.D2",
    border_gp = gpar(lwd = 1),
    show_heatmap_legend = TRUE,
    row_title_rot = 0
  )

  #draw the heatmap with the bottom annotation to file
  CairoPDF(outpath, width = plot_width, height = plot_height)
  draw(hm_signature)
  dev.off()
}
  