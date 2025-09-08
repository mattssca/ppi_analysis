#source hm order fun
source("functions/get_hm_order.R")

nodes_object = erbb_nodes_small

#extended gene network (for plotting)
genes_to_plot = nodes_object$expanded_genes

#get subtype colors
col = list(`Subtype` = lund_colors$lund_colors)

#get labels
pred_lab = subtype_vector_class_7_stripped

#get split
split = factor(pred_lab, levels = c("UroA", "UroB", "UroC", "GU", "BaSq"))

#heatmap colors
col_fun = circlize::colorRamp2(c(-2, 0, 2), c("#4DF76F", "black", "#F74D4D"))

#create signature score
#subset the expression matrix to only the genes in genes_to_plot
expr_subset <- expr_z_sub[rownames(expr_z_sub) %in% genes_to_plot, , drop = FALSE]
expr_z_sub = as.matrix(expr_z_sub)

#calculate the mean for each sample across the selected genes
signature_score <- apply(expr_subset, 2, median, na.rm = TRUE)

#signature color
sig_min <- min(signature_score)
sig_max <- max(signature_score)
col_signature = colorRamp2(c(sig_min, sig_max), c("#1E93AB", "#E62727"))

#get degree as an annotation track
degree_vec <- igraph::degree(nodes_object$g)
degree_vec <- degree_vec[genes_to_plot]

deg_min <- min(degree_vec, na.rm = TRUE)
deg_max <- max(degree_vec, na.rm = TRUE)
col_degree <- colorRamp2(c(deg_min, deg_max), c("#BCA88D", "#3E3F29"))

#bind the color objects into a list
col_list = list(Subtype = lund_colors$lund_colors,
           `Signature Score` = col_signature,
           Degree = col_degree)

#get sample order
sample_order = get_sample_order(these_predictions = sub_predicted)


#draw top annotation
hm_annotation_top = HeatmapAnnotation(`Subtype` = pred_lab,
                                      annotation_name_side = "left",
                                      gap = unit(1, "mm"),
                                      simple_anno_size = unit(4, "mm"),
                                      simple_anno_size_adjust = TRUE,
                                      col = col_list,
                                      show_legend = FALSE,
                                      border = TRUE,
                                      annotation_name_gp = gpar(fontsize = 8))

#draw bottom annotation
hm_annotation_bottom = HeatmapAnnotation(`Signature Score` = signature_score,
                                         annotation_name_side = "left",
                                         gap = unit(1, "mm"),
                                         simple_anno_size = unit(4, "mm"),
                                         simple_anno_size_adjust = TRUE,
                                         col = col_list,
                                         show_legend = TRUE,
                                         border = TRUE,
                                         annotation_name_gp = gpar(fontsize = 8))

hm_annotation_left <- rowAnnotation(Degree = degree_vec,
                                    col = col_list,
                                    show_annotation_name = FALSE)

#draw signature heatmap
hm_signature = Heatmap(expr_z_sub[genes_to_plot,,drop = FALSE],
                                             name = "Expression",
                                             top_annotation = hm_annotation_top, 
                                             bottom_annotation = hm_annotation_bottom,
                                             left_annotation = hm_annotation_left,
                                             col = col_fun,
                                             column_split = split,
                                             column_title = "ERBB Extended PPI Network",
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
                                             row_names_gp = gpar(fontsize = 8,
                                                                 fontface = ifelse(rownames(expr_z_sub[genes_to_plot,,drop = FALSE]) %in% nodes_object$signature_genes, "bold", "plain"),
                                                                 col = ifelse(rownames(expr_z_sub[genes_to_plot,,drop = FALSE]) %in% nodes_object$signature_genes, "#E62727", "black")),
                                             row_title_gp = gpar(fontsize = 10),
                                             cluster_rows = TRUE,
                                             clustering_distance_rows = "spearman",
                                             clustering_method_rows = "ward.D2",
                                             border_gp = gpar(lwd = 1),
                                             show_heatmap_legend = TRUE,
                                             row_title_rot = 0)


#export plot
CairoPDF("viz/tat1/ERBB/complete_analysis/figures/small/score_high/extended_signature_heatmap.pdf", width = 20, height = 20)
draw(hm_signature)
dev.off()




