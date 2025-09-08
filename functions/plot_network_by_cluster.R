#' Plot PPI Network Colored by Cluster Membership
#'
#' Visualizes a PPI network using ggraph, coloring nodes by cluster assignment.
#'
#' @param nodes_object List. Output from expand_and_plot_signature_network (must contain $g).
#' @param clusters Named vector. Cluster assignments for each gene (names = gene symbols, values = cluster numbers).
#' @param cluster_colors Named vector. Colors for each cluster (default: 3 clusters).
#' @param layout_coords Matrix. Node layout coordinates from expand_and_plot_signature_network (optional).
#' @param plot_title Character. Title for the plot.
#'
#' @return A ggplot object (network plot).
#' @export
#' 
plot_network_by_cluster <- function(nodes_object, 
                                    clusters, 
                                    cluster_colors = c("1" = "#77BEF0", "2" = "#FFCB61", "3" = "#EA5B6F"), 
                                    node_size = 17,
                                    pdf_width = 15,
                                    pdf_height = 10,
                                    output_path, 
                                    plot_title = "Network Colored by Cluster Membership"){
  library(igraph)
  library(ggraph)
  library(tidygraph)
  library(ggplot2)
  
  
  if(length(unique(clusters)) > length(unique(cluster_colors))) {
    stop(
      paste0(
        "Not enough colors provided for annotating clusters: ",
        length(unique(clusters)), " clusters detected, but only ",
        length(unique(cluster_colors)), " colors provided."
      )
    )
  }
  
  #convert igraph object to tidygraph
  tg <- as_tbl_graph(nodes_object$g)
  this_layout = nodes_object$layout_coords

  #add cluster membership as a node attribute
  tg <- tg %>% mutate(cluster = clusters[name])
  
  p <- ggraph(tg, layout = this_layout) +
    geom_edge_link(color = "#969696", alpha = 1, width = 0.3) +
    geom_node_point(aes(color = as.factor(cluster)), size = node_size, show.legend = FALSE) +
    geom_node_text(aes(label = name), repel = TRUE, size = 3) +
    scale_color_manual(values = cluster_colors, name = "Cluster") +
    geom_label(aes(x = x, y = y, label = name), fill = "white", color = "black", size = 6, fontface = "bold", label.size = 0.2) +
  theme_bw(base_size = 15) +
    theme(
      axis.title.x = element_blank(), 
      axis.title.y = element_blank(),
      axis.text.x = element_blank(), 
      panel.grid = element_blank(),
      axis.text.y = element_blank(), 
      axis.ticks = element_blank(),
      plot.title.position = "plot",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    ) +
    labs(title = plot_title)
  
    CairoPDF(output_path, width = pdf_width, height = pdf_height)
    print(p)
    dev.off()
}
