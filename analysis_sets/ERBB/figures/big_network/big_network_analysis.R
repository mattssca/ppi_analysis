#load data
load("analysis_sets/uroscanseq_data.Rdata")

#load packages
library(ggplot2)
library(dplyr)
library(igraph)
library(STRINGdb)
library(ggraph)
library(Cairo)
library(biomaRt)
library(ComplexHeatmap)
library(circlize)
library(patchwork)

#source functions
source("functions/expand_and_plot_signature_network.R")
source("functions/plot_network_by_cluster.R")
source("functions/plot_node_heatmap.R") 
source("functions/get_hm_order.R")
source("functions/plot_signature_network_heatmap.R")
source("functions/run_anova_ppi.R")

#get lundtax objects
lund_colors = LundTax2023Classifier::lund_colors$lund_colors
signatures = LundTax2023Classifier::signatures

#create nodes object
erbb_nodes = expand_and_plot_signature_network(expr_data = uroscanseq_data$expr_df,
                                               return_data = TRUE,
                                               expr_summary = "mean",
                                               max_added_genes = 300, 
                                               layout_method = "kk",
                                               min_degree = 1, 
                                               string_score_threshold = 500,  
                                               subtype_vector = uroscanseq_data$subtype_7_vector, 
                                               signature_list = uroscanseq_data$signature_genes, 
                                               signature_name = "ERBB")

#plot ppi network
expand_and_plot_signature_network(expr_data = uroscanseq_data$expr_df,
                                  node_size = 10,
                                  plot_width = 30,
                                  plot_height = 30, 
                                  color_scale_global = FALSE,
                                  return_data = FALSE,
                                  expr_summary = "mean",
                                  node_color = "lund",
                                  node_degree = TRUE,
                                  theme = "light",
                                  show_labels = TRUE,
                                  max_added_genes = 300, 
                                  layout_method = "kk",
                                  min_degree = 1, 
                                  string_score_threshold = 500,  
                                  subtype_vector = uroscanseq_data$subtype_7_vector, 
                                  signature_list = uroscanseq_data$signature_genes, 
                                  signature_name = "ERBB",
                                  verbose = TRUE, 
                                  out_dir = "analysis_sets/ERBB/figures/big_network/")

#plot node heatmap
erbb_cluster = plot_node_heatmap(nodes_object = erbb_nodes, 
                                 output_path = "analysis_sets/ERBB/figures/big_network/ERBB_nodes_heatmap.pdf", 
                                 plot_title = "ERBB PPI Heatmap", 
                                 top_hubs = 50,
                                 cell_width = 2, 
                                 pdf_width = 30, 
                                 pdf_height = 10, 
                                 column_split = 3)

#plot network by cluster
plot_network_by_cluster(nodes_object = erbb_nodes,
                        node_size = 20,
                        cluster_colors = c("1" = "#77BEF0", "2" = "#FFCB61", "3" = "#EA5B6F"),
                        pdf_width = 30, 
                        pdf_height = 30, 
                        clusters = erbb_cluster, 
                        plot_title = "ERBB PPI Network by Clusters", 
                        output_path = "analysis_sets/ERBB/figures/big_network/ERBB_PPI_clusters.pdf")

#plot node HUB heatmap
erbb_cluster_hubs = plot_node_heatmap(nodes_object = erbb_nodes,
                                      only_hubs = TRUE, 
                                      output_path = "analysis_sets/ERBB/figures/big_network/ERBB_nodes_HUB_heatmap.pdf", 
                                      plot_title = "ERBB PPI Heatmap - Top 10 HUBS", 
                                      top_hubs = 50,
                                      pdf_width = 20, 
                                      cell_width = 5,
                                      column_split = 3)

#create extended signature heatmap (per sample)
plot_signature_network_heatmap(nodes_object = erbb_nodes, 
                               expr_data = uroscanseq_data$expr_df, 
                               subtype_vector = uroscanseq_data$subtype_7_vector,
                               lund_colors = lund_colors,
                               only_hubs = FALSE, 
                               global_signature_cols = FALSE,
                               top_hubs = 50,
                               outpath = "analysis_sets/ERBB/figures/big_network/ERBB_signature_heatmap.pdf", 
                               plot_width = 70,
                               plot_height = 50,
                               cluster_object = erbb_cluster,
                               plot_title = "ERBB Signature (Extended) Heatmap")

 plot_signature_network_heatmap(nodes_object = erbb_nodes, 
                               expr_data = uroscanseq_data$expr_df, 
                               subtype_vector = uroscanseq_data$subtype_7_vector,
                               lund_colors = lund_colors,
                               only_hubs = TRUE, 
                               top_hubs = 10,
                               outpath = "analysis_sets/ERBB/figures/big_network/ERBB_signature_heatmap_HUBS_only.pdf", 
                               plot_width = 30,
                               plot_height = 30,
                               cluster_object = erbb_cluster_hubs,
                               plot_title = "ERBB Signature (Top 10 HUBS) Heatmap")

#run ANOVA
erbb_results = run_anova_ppi(expr_data = uroscanseq_data$expr_df, 
                             extended_genes = erbb_nodes$expanded_genes, 
                             subtype_vector = uroscanseq_data$subtype_7_vector,
                             node_metrics = erbb_nodes$node_metrics, 
                             subtype_class = "7_class", 
                             sig_p = 0.05,
                             min_diff = 1)

CairoPDF("analysis_sets/ERBB/figures/big_network/ERBB_anova_volcano.pdf", width = 15, height = 10)
print(erbb_results$volcano_plot)
dev.off()

CairoPDF("analysis_sets/ERBB/figures/big_network/ERBB_facet_expression_volcano.pdf", width = 20, height = 20)
print(erbb_results$facet_plot)
dev.off()

#redraw heatmap with blacklsited genes
plot_blacklist_network_heatmap(nodes_object = erbb_nodes, 
                               expr_data = uroscanseq_data$expr_df, 
                               subtype_vector = uroscanseq_data$subtype_7_vector,
                               lund_colors = lund_colors,
                               blacklist_genes = erbb_results$blacklist_genes,
                               cut_tree = 7,
                               subtype_class = "7_class",
                               outpath = "analysis_sets/ERBB/figures/big_network/ERBB_blacklist_signature_heatmap.pdf", 
                               plot_width = 20,
                               plot_height = 30,
                               plot_title = "ERBB Signature Blacklisted Genes Removed")
 
