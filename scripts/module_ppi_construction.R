#load packages
library(igraph)
library(ggraph)
library(dplyr)

#subset to relevant genes
erbb_genes = as.character(erbb_big_network$go_enrichment$gene)

#detect communities (e.g., Louvain)
comm <- cluster_louvain(erbb_big_network$igraph_object)
modules <- membership(comm)

#filter modules to only include genes in erbb_genes
modules_in_erbb <- modules[names(modules) %in% erbb_genes]

#split genes by module
genes_by_module <- split(names(modules_in_erbb), modules_in_erbb)

#create a PPI Network for Each Module
ppi_networks_by_module <- lapply(genes_by_module, function(genes){
  induced_subgraph(erbb_big_network$igraph_object, vids = genes)
})


#calculate mean expression for each subtype
mean_expr_by_subtype <- sapply(unique(erbb_big_network$subtype_vector), function(st) {
  samples <- names(erbb_big_network$subtype_vector)[erbb_big_network$subtype_vector == st]
  rowMeans(erbb_big_network$expr_data[, samples, drop = FALSE], na.rm = TRUE)
})


#plot each module and subtype
for (mod in names(genes_by_module)) {
  genes <- genes_by_module[[mod]]
  g <- induced_subgraph(erbb_big_network$igraph_object, vids = genes)
  
  #precompute layout
  layout_mat <- layout_with_fr(g, niter = 1000)
  layout_df <- as.data.frame(layout_mat)
  colnames(layout_df) <- c("x", "y")
  layout_df$name <- V(g)$name
  
  for (st in colnames(mean_expr_by_subtype)) {
    expr_vals <- mean_expr_by_subtype[genes, st]
    V(g)$mean_expr <- expr_vals[V(g)$name]
    subtype_range <- range(expr_vals, na.rm = TRUE)
    
    # Use the precomputed layout
    p <- ggraph(g, layout = "manual", x = layout_df$x, y = layout_df$y) +
      geom_edge_link(color = "black", alpha = 1, width = 0.3) +
      geom_node_point(aes(color = mean_expr), size = 20, show.legend = FALSE) +
      scale_color_gradientn(
        colors = c("#4DF76F", "black", "#F74D4D"),
        values = scales::rescale(c(subtype_range[1], 0, subtype_range[2])),
        limits = subtype_range,
        na.value = "grey80",
        name = "Mean Expr"
      ) +
      geom_label(aes(x = x, y = y, label = name), fill = "white", color = "black", size = 4, fontface = "bold", label.size = 0.2) +
      theme_bw(base_size = 15) +
      theme(
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(), 
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        plot.title.position = "plot",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5)
      ) +
      labs(title = paste("Module", mod, "-", st))
    
    ggsave(filename = paste0("analysis_sets/ERBB/figures/big_network/module_ppi/module_", mod, "_", st, ".png"), plot = p, width = 10, height = 7, dpi = 600)
  }
}
