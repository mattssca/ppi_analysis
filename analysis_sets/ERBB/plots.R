#source scripts
source("functions/plot_network_by_cluster.R")
source("functions/plot_node_heatmap.R")
source("analysis/expand_and_plot_signature_network.R")


############# PLOT PPI NETOWRK ############# 
expand_and_plot_signature_network(expr_data = expr_z_sub,
                                  node_size = 20,
                                  plot_width = 8,
                                  plot_height = 8, 
                                  color_scale_global = FALSE,
                                  return_data = FALSE,
                                  expr_summary = "mean",
                                  node_color = "lund",
                                  node_degree = TRUE,
                                  theme = "light",
                                  show_labels = TRUE,
                                  max_added_genes = 20, 
                                  layout_method = "kk",
                                  min_degree = 1, 
                                  string_score_threshold = 700,  
                                  subtype_vector = subtype_vector_class_7_stripped, 
                                  signature_list = lundtax_signatures, 
                                  signature_name = "ERBB",
                                  verbose = TRUE, 
                                  out_dir = "viz/tat1/ERBB/complete_analysis/figures/small/score_high/")

 expand_and_plot_signature_network(expr_data = expr_z_sub,
                                  node_size = 25,
                                  plot_width = 20,
                                  plot_height = 15, 
                                  color_scale_global = FALSE,
                                  return_data = FALSE,
                                  expr_summary = "mean",
                                  node_color = "lund",
                                  node_degree = FALSE,
                                  theme = "light",
                                  show_labels = TRUE,
                                  max_added_genes = 200, 
                                  layout_method = "kk",
                                  min_degree = 1, 
                                  string_score_threshold = 500,  
                                  subtype_vector = subtype_vector_class_7_stripped, 
                                  signature_list = lundtax_signatures, 
                                  signature_name = "ERBB",
                                  verbose = TRUE, 
                                  out_dir = "viz/tat1/ERBB/complete_analysis/figures/large/")

############# PLOT NODE HEATMAP #############
erbb_cluster_small = plot_node_heatmap(nodes_object = erbb_nodes_small, 
                                       output_path = "viz/tat1/ERBB/complete_analysis/figures/small/score_high/node_heatmap.pdf", 
                                       plot_title = "ERBB Small", 
                                       top_hubs = 10, 
                                       column_split = 3)


#remove NA genes
erbb_nodes_large$node_metrics <- erbb_nodes_large$node_metrics[!erbb_nodes_large$node_metrics$name %in% c("IFNA1", "SOX2"), ]

erbb_cluster_large = plot_node_heatmap(nodes_object = erbb_nodes_large, 
                                       output_path = "viz/tat1/ERBB/complete_analysis/figures/large/node_heatmap.pdf", 
                                       plot_title = "ERBB Large", 
                                       top_hubs = 10,
                                       pdf_width = 35, 
                                       cell_width = 4,
                                       column_split = 3)

########## PLOT NETWORK BY CLUSTER ########## 
plot_network_by_cluster(nodes_object = erbb_nodes_small,
                        node_size = 20,
                        pdf_width = 8, 
                        pdf_height = 8, 
                        clusters = erbb_cluster_small, 
                        plot_title = "ERBB Small", 
                        output_path = "viz/tat1/ERBB/complete_analysis/figures/small/score_high/erbb_clsuters.pdf")

plot_network_by_cluster(nodes_object = erbb_nodes_large, 
                        clusters = erbb_cluster_large,
                        node_size = 25,
                        pdf_width = 20, 
                        pdf_height = 15, 
                        plot_title = "ERBB Large", 
                        output_path = "viz/tat1/ERBB/complete_analysis/figures/large/erbb_clsuters.pdf")

############# PLOT NODE HEATMAP ONLY HUBS #############
erbb_cluster_small_hubs = plot_node_heatmap(nodes_object = erbb_nodes_small,
                                            only_hubs = TRUE, 
                                            output_path = "viz/tat1/ERBB/complete_analysis/figures/small/score_high/node_heatmap_hubs.pdf", 
                                            plot_title = "ERBB Small Hubs", 
                                            top_hubs = 10,
                                            pdf_width = 10, 
                                            cell_width = 10,
                                            column_split = 2)

erbb_cluster_large_hubs = plot_node_heatmap(nodes_object = erbb_nodes_large, 
                                            only_hubs = TRUE,
                                            output_path = "viz/tat1/ERBB/complete_analysis/figures/large/node_heatmap_hubs.pdf", 
                                            plot_title = "ERBB Large Hubs", 
                                            top_hubs = 10,
                                            pdf_width = 10, 
                                            cell_width = 10,
                                            column_split = 2)
#LundTax heatmap
#subset the 'data' data frame
sub_predicted = lundtax_predict_sub(this_data = expr_z_sub, 
                                    adjust = TRUE, 
                                    impute = TRUE, 
                                    include_data = TRUE, 
                                    log_transform = FALSE)

plot_hm_signatures(these_predictions = sub_predicted, 
                   norm = FALSE, 
                   ann_height = 5,
                   subtype_annotation = "7_class",
                   out_path = "viz/tat1/ERBB/complete_analysis/figures/", 
                   out_format = "pdf", 
                   plot_scores = TRUE, 
                   plot_title = "ERBB Signature")

# Convert to data frame
subtype_df <- data.frame(
  subtype = as.character(subtype_vector_class_7_stripped)
)

# Set factor levels to desired order
subtype_df$subtype <- factor(subtype_df$subtype, levels = c("UroA", "UroB", "UroC", "GU", "BaSq"))

# Define your custom colors
subtype_colors <- c(
  "UroA" = "#3cb44b",
  "UroB" = "#8B1A1A",
  "UroC" = "#006400",
  "GU"   = "#4363d8",
  "BaSq" = "#CD2626"
)

# Count samples per subtype for annotation
subtype_counts <- as.data.frame(table(subtype_df$subtype))
colnames(subtype_counts) <- c("subtype", "count")

# Plot with counts above bars
ggplot(subtype_df, aes(x = subtype, fill = subtype)) +
  geom_bar() +
  geom_text(data = subtype_counts, aes(x = subtype, y = count, label = count), 
            vjust = -0.5, color = "black", size = 3) +
  scale_fill_manual(values = subtype_colors) +
  labs(title = "Number of Samples per Subtype",
       x = "Subtype",
       y = "Number of Samples") +
  theme_minimal()
