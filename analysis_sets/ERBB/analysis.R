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

#source functions
source("functions/expand_and_plot_signature_network.R")
source("functions/plot_network_by_cluster.R")
source("functions/plot_node_heatmap.R") 
source("functions/get_hm_order.R")
source("functions/plot_signature_network_heatmap.R")

#get lundtax objects
lund_colors = LundTax2023Classifier::lund_colors$lund_colors
signatures = LundTax2023Classifier::signatures

#create nodes object
erbb_nodes = expand_and_plot_signature_network(expr_data = uroscanseq_data$expr_df,
                                               return_data = TRUE,
                                               expr_summary = "mean",
                                               max_added_genes = 30, 
                                               layout_method = "kk",
                                               min_degree = 1, 
                                               string_score_threshold = 700,  
                                               subtype_vector = uroscanseq_data$subtype_7_vector, 
                                               signature_list = uroscanseq_data$signature_genes, 
                                               signature_name = "ERBB")

#plot ppi netowrk
expand_and_plot_signature_network(expr_data = uroscanseq_data$expr_df,
                                  node_size = 20,
                                  plot_width = 10,
                                  plot_height = 10, 
                                  color_scale_global = FALSE,
                                  return_data = FALSE,
                                  expr_summary = "mean",
                                  node_color = "lund",
                                  node_degree = TRUE,
                                  theme = "light",
                                  show_labels = TRUE,
                                  max_added_genes = 30, 
                                  layout_method = "kk",
                                  min_degree = 1, 
                                  string_score_threshold = 700,  
                                  subtype_vector = uroscanseq_data$subtype_7_vector, 
                                  signature_list = uroscanseq_data$signature_genes, 
                                  signature_name = "ERBB",
                                  verbose = TRUE, 
                                  out_dir = "analysis_sets/ERBB/figures/")

#plot node heatmap
erbb_cluster = plot_node_heatmap(nodes_object = erbb_nodes, 
                                 output_path = "analysis_sets/ERBB/figures/ERBB_nodes_heatmap.pdf", 
                                 plot_title = "ERBB PPI Heatmap", 
                                 top_hubs = 10, 
                                 column_split = 3)

#plot network by cluster
plot_network_by_cluster(nodes_object = erbb_nodes,
                        node_size = 20,
                        cluster_colors = c("1" = "#77BEF0", "2" = "#FFCB61", "3" = "#EA5B6F"),
                        pdf_width = 10, 
                        pdf_height = 10, 
                        clusters = erbb_cluster, 
                        plot_title = "ERBB PPI Network by Clusters", 
                        output_path = "analysis_sets/ERBB/figures/ERBB_PPI_clusters.pdf")

#plot node HUB heatmap
erbb_cluster_hubs = plot_node_heatmap(nodes_object = erbb_nodes,
                                      only_hubs = TRUE, 
                                      output_path = "analysis_sets/ERBB/figures/ERBB_nodes_HUB_heatmap.pdf", 
                                      plot_title = "ERBB PPI Heatmap - Top 10 HUBS", 
                                      top_hubs = 10,
                                      pdf_width = 10, 
                                      cell_width = 10,
                                      column_split = 2)

#create extended signature heatmap (per sample)
plot_signature_network_heatmap(nodes_object = erbb_nodes, 
                               expr_data = uroscanseq_data$expr_df, 
                               subtype_vector = uroscanseq_data$subtype_7_vector,
                               lund_colors = lund_colors,
                               only_hubs = FALSE, 
                               global_signature_cols = FALSE,
                               top_hubs = 10,
                               outpath = "analysis_sets/ERBB/figures/ERBB_signature_heatmap.pdf", 
                               plot_width = 20,
                               plot_height = 10,
                               cluster_object = erbb_cluster,
                               plot_title = "ERBB Signature (Extended) Heatmap")

plot_signature_network_heatmap(nodes_object = erbb_nodes, 
                               expr_data = uroscanseq_data$expr_df, 
                               subtype_vector = uroscanseq_data$subtype_7_vector,
                               lund_colors = lund_colors,
                               only_hubs = TRUE, 
                               top_hubs = 10,
                               outpath = "analysis_sets/ERBB/figures/ERBB_signature_heatmap_HUBS_only.pdf", 
                               plot_width = 20,
                               plot_height = 10,
                               cluster_object = erbb_cluster_hubs,
                               plot_title = "ERBB Signature (Top 10 HUBS) Heatmap")
 
