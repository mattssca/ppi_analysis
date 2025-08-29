#' Plot Extended PPI Network Heatmap
#'
#' Generates a ComplexHeatmap of the extended PPI network with annotation tracks for group, degree, 
#' and community.
#'
#' @param nodes_object List. Output from expand_and_plot_signature_network (must contain node_metrics, signature_genes, etc.)
#' @param output_path Character. Path to save the PDF output.
#' @param plot_title Character. Title for the heatmap plot.
#' @param column_split Integer. Number of clusters to cut the column dendrogram (default: 3).
#' @param cell_width Numeric. Width of each heatmap cell in millimeters (default: 10).
#' @param cell_height Numeric. Height of each heatmap cell in millimeters (default: 10).
#' @param pdf_width Numeric. Width of the output PDF in inches (default: 20).
#' @param pdf_height Numeric. Height of the output PDF in inches (default: 10).
#' @param top_hubs Integer. Number of top hub genes to highlight in the annotation track (default: 5).
#' 
#' @import ComplexHeatmap Cairo circlize viridis
#'
#' @return No return value. The heatmap is saved to the specified output_path.
#' @export
#' 
plot_node_heatmap <- function(nodes_object, 
                              output_path, 
                              plot_title = NULL,
                              column_split = 3,
                              cell_width = NULL,
                              cell_height = NULL,
                              pdf_width = 20,
                              pdf_height = 10,
                              top_hubs = 5){
  library(ComplexHeatmap)
  library(Cairo)
  library(circlize)
  library(viridis)
  
  node_metrics <- nodes_object$node_metrics
  lund_genes <- LundTax2023Classifier::gene_list$hgnc_symbol
  gene_names <- nodes_object$expanded_genes
  
  # Order subtypes as UroA, UroB, UroC, GU, BaSq
  desired_order <- c("mean_expr_UroA", "mean_expr_UroB", "mean_expr_UroC", "mean_expr_GU", "mean_expr_BaSq")
  expr_cols <- grep("^mean_expr_", colnames(node_metrics), value = TRUE)
  expr_cols <- intersect(desired_order, expr_cols)
  expr_matrix <- as.matrix(node_metrics[, expr_cols])
  rownames(expr_matrix) <- node_metrics$name
  colnames(expr_matrix) <- gsub("^mean_expr_", "", colnames(expr_matrix))
  
  # Transpose to flip heatmap (genes as columns)
  expr_matrix_t <- t(scale(expr_matrix))
  
  # Create annotation for columns (now genes)
  gene_group <- rep("Neighbors (Extended)", ncol(expr_matrix_t))
  colnames(expr_matrix_t) <- rownames(expr_matrix) # ensure column names match gene names
  gene_names <- colnames(expr_matrix_t)
  gene_group[gene_names %in% lund_genes] <- "Neighbors (Lund)"
  gene_group[gene_names %in% nodes_object$signature_genes] <- "Lund Signature Genes"
  
  group_colors <- c(
    "Neighbors (Extended)" = "#F7AD45",
    "Neighbors (Lund)" = "#BB3E00",
    "Lund Signature Genes" = "#657C6A"
  )
  
  degree_colors <- colorRamp2(
    c(min(node_metrics$degree), max(node_metrics$degree)),
    c("#BCA88D", "#3E3F29")
  )
  
  
  # Identify top N hubs by degree
  top_hub_names <- node_metrics$name[order(node_metrics$degree, decreasing = TRUE)][1:top_hubs]
  hub_status <- ifelse(gene_names %in% top_hub_names, "Top Hub", "Other")
  hub_colors <- c("Top Hub" = "#E62727", "Other" = "#1E93AB")
  
  col_ha <- HeatmapAnnotation(
    annotation_name_side = "left", 
    `Gene Group` = gene_group,
    Degree = node_metrics$degree[match(gene_names, node_metrics$name)],
    Hub = hub_status,
    col = list(
      `Gene Group` = group_colors,
      Degree = degree_colors,
      Hub = hub_colors
    ),
    show_annotation_name = TRUE
  )
  
  heatmap_colors <- colorRampPalette(c("#4DF76F", "black", "#F74D4D"))(100)
  
  #draw heatmap
  # Calculate cell width and height so that total width = 350mm and height = 50mm by default
  n_cols <- ncol(expr_matrix_t)
  n_rows <- nrow(expr_matrix_t)
  if (is.null(cell_width)) cell_width <- 350 / n_cols
  if (is.null(cell_height)) cell_height <- 50 / n_rows
  
  ht = Heatmap(
    expr_matrix_t, 
    name = "Expression",
    col = heatmap_colors,
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    column_split = column_split,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 14), 
    row_names_side = "left", 
    column_names_gp = gpar(fontsize = 12),      
    column_names_side = "bottom",              
    rect_gp = gpar(col = NA),                   
    top_annotation = col_ha,
    show_heatmap_legend = FALSE, 
    width = unit(cell_width * n_cols, "mm"),
    height = unit(cell_height * n_rows, "mm"),
    column_title = plot_title
  )
  
  CairoPDF(output_path, width = pdf_width, height = pdf_height)
  draw(ht)
  dev.off()
}
