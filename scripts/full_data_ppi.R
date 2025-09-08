# Full Cell Cycle PPI Network with Subtype-Specific Expression Coloring
# This script creates a single PPI network using all cell cycle genes
# and colors nodes based on which subtype they have highest expression in

# Load required libraries
library(igraph)
library(ggraph)
library(tidygraph)
library(ggplot2)
library(dplyr)

# Load data
load("data/expr/expr_cell_cycle.Rdata")  # Full cell cycle expression matrix
load("data/expr/expr_subtypes.Rdata")    # Subtype-specific expression data

# Source the PPI function
source("functions/get_ppi.R")

# Define Lund colors
lund_colors <- c(
  "uro" = "#3cb44b",
  "gu" = "#4363d8", 
  "basq" = "#CD2626",
  "mes" = "#f58231",
  "scne" = "#A020F0"
)

# Step 1: Calculate mean expression for each gene in each subtype
gene_subtype_expression <- data.frame()

for(subtype in names(expr_subtypes)) {
  # Calculate mean expression for each gene in this subtype
  mean_expr <- rowMeans(expr_subtypes[[subtype]], na.rm = TRUE)
  
  temp_df <- data.frame(
    gene = names(mean_expr),
    subtype = subtype,
    mean_expression = mean_expr
  )
  
  gene_subtype_expression <- rbind(gene_subtype_expression, temp_df)
}

# Step 2: For each gene, find which subtype has the highest expression
gene_dominant_subtype <- gene_subtype_expression %>%
  group_by(gene) %>%
  slice_max(mean_expression, n = 1) %>%
  ungroup() %>%
  select(gene, dominant_subtype = subtype, max_expression = mean_expression)

# Check the results
cat("Genes per dominant subtype:\n")
table(gene_dominant_subtype$dominant_subtype)

# Step 3: Create PPI network using the full cell cycle expression data
cat("\nCreating PPI network from full cell cycle expression data...\n")
full_ppi_network <- run_stringdb(expr_data = expr_cell_cycle)

cat("Network created with", vcount(full_ppi_network), "nodes and", ecount(full_ppi_network), "edges\n")

# Step 4: Add subtype dominance information to the network
# Get node names from the network
network_genes <- V(full_ppi_network)$name

# Match with our dominant subtype data
gene_colors <- gene_dominant_subtype[match(network_genes, gene_dominant_subtype$gene), ]

# Add attributes to the network
V(full_ppi_network)$dominant_subtype <- gene_colors$dominant_subtype
V(full_ppi_network)$max_expression <- gene_colors$max_expression
V(full_ppi_network)$node_color <- lund_colors[gene_colors$dominant_subtype]

# Handle any missing values (genes not found in subtype data)
missing_genes <- is.na(V(full_ppi_network)$dominant_subtype)
if(any(missing_genes)) {
  V(full_ppi_network)$dominant_subtype[missing_genes] <- "unknown"
  V(full_ppi_network)$node_color[missing_genes] <- "gray"
  cat("Warning:", sum(missing_genes), "genes not found in subtype data, colored gray\n")
}

# Step 5: Create the visualization function
create_full_ppi_plot <- function(top_n = 100, layout_type = "fr", show_labels = TRUE, highlight_top_hubs = 10) {
  
  # Get top N nodes by degree for cleaner visualization
  degrees <- degree(full_ppi_network)
  top_nodes <- names(sort(degrees, decreasing = TRUE)[1:top_n])
  subgraph <- induced_subgraph(full_ppi_network, top_nodes)
  
  # Identify top hub genes in the subgraph
  subgraph_degrees <- degree(subgraph)
  top_hub_genes <- names(sort(subgraph_degrees, decreasing = TRUE)[1:highlight_top_hubs])
  
  # Convert to tidygraph
  tidy_subgraph <- as_tbl_graph(subgraph) %>%
    activate(nodes) %>%
    mutate(
      degree = degree(subgraph),
      # Get the attributes we added
      dominant_subtype = V(subgraph)$dominant_subtype,
      max_expression = V(subgraph)$max_expression,
      node_color = V(subgraph)$node_color,
      # Add hub highlighting
      is_top_hub = name %in% top_hub_genes
    )
  
  # Create the plot
  p <- ggraph(tidy_subgraph, layout = layout_type) +
    geom_edge_link(alpha = 0.3, color = "gray70", width = 0.5) +
    geom_node_point(aes(size = degree, alpha = max_expression, 
                        stroke = ifelse(is_top_hub, 2, 0.5)), 
                    color = "black", fill = V(subgraph)$node_color, 
                    shape = 21) +
    scale_alpha_continuous(range = c(0.4, 1.0), name = "Max Expression") +
    scale_size_continuous(range = c(3, 12), name = "Degree") +
    labs(title = "Cell Cycle PPI Network - Colored by Dominant Subtype",
         subtitle = paste("Top", top_n, "nodes by degree. Node color = subtype with highest expression. Thick outline =", highlight_top_hubs, "top hubs"),
         caption = "Green=Uro, Blue=GU, Red=BaSq, Orange=Mes, Purple=ScNE") +
    theme_graph() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, family = "Arial"),
      plot.subtitle = element_text(size = 12, hjust = 0.5, family = "Arial"),
      plot.caption = element_text(size = 10, hjust = 0.5, family = "Arial"),
      legend.position = "none"
    )
  
  # Add labels if requested
  if (show_labels) {
    p <- p + geom_node_text(aes(label = name), size = 2.5, repel = TRUE, 
                           max.overlaps = 30, color = "black", fontface = "bold",
                           family = "Arial")
  }
  
  return(p)
}

# Step 6: Create and display the plots
full_plot_clean <- create_full_ppi_plot(top_n = 50, 
                                        show_labels = TRUE, 
                                        highlight_top_hubs = 10,
                                        layout_type = "fr")
print(full_plot_clean)


# Save the plot with Cairo PDF at 20x15 inches
cairo_pdf("full_cell_cycle_ppi_network_50.pdf", 
          width = 10, height = 7, 
          pointsize = 12)
print(full_plot_clean)
dev.off()

# Step 7: Analysis summary
cat("\n=== ANALYSIS SUMMARY ===\n")
cat("Total genes in network:", vcount(full_ppi_network), "\n")
cat("Total edges in network:", ecount(full_ppi_network), "\n")
cat("\nGenes by dominant subtype:\n")
subtype_summary <- table(V(full_ppi_network)$dominant_subtype)
print(subtype_summary)

# Show top hub genes for each dominant subtype
cat("\n=== TOP HUB GENES BY DOMINANT SUBTYPE ===\n")
network_analysis <- data.frame(
  gene = V(full_ppi_network)$name,
  degree = degree(full_ppi_network),
  dominant_subtype = V(full_ppi_network)$dominant_subtype,
  max_expression = V(full_ppi_network)$max_expression
) %>%
  arrange(desc(degree))

for(subtype in names(lund_colors)) {
  subtype_genes <- network_analysis %>%
    filter(dominant_subtype == subtype) %>%
    slice_head(n = 5)
  
  if(nrow(subtype_genes) > 0) {
    cat("\nTop", subtype, "genes (", lund_colors[subtype], "):\n")
    print(subtype_genes)
  }
}
