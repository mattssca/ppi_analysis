library(pheatmap)
expr_cols <- grep("^mean_expr_", colnames(node_metrics), value = TRUE)
expr_matrix <- as.matrix(node_metrics[, expr_cols])
rownames(expr_matrix) <- node_metrics$name
pheatmap(expr_matrix, cluster_rows = TRUE, cluster_cols = TRUE, scale = "row")


library(reshape2)
df <- melt(node_metrics[, c("name", expr_cols)], id.vars = "name")
library(ggplot2)
ggplot(df, aes(x = variable, y = value, fill = name)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal()


ggplot(node_metrics, aes(x = degree, y = mean_expr_GU)) +
  geom_point() +
  theme_minimal() +
  labs(x = "Degree", y = "Mean Expression (GU)")


library(ggplot2)
# Example: Count genes per community and subtype
bubble_data <- reshape2::melt(table(node_metrics$community, apply(expr_matrix, 1, function(x) names(which.max(x)))))
colnames(bubble_data) <- c("community", "subtype", "count")
ggplot(bubble_data, aes(x = community, y = subtype, size = count)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(title = "Subtype Enrichment in Communities")

library(UpSetR)
# Create a list of genes for each subtype
gene_lists <- lapply(expr_cols, function(col) node_metrics$name[node_metrics[[col]] > 1])
names(gene_lists) <- sub("mean_expr_", "", expr_cols)
upset(fromList(gene_lists), order.by = "freq")


library(ggplot2)
dot_data <- reshape2::melt(node_metrics[, c("name", expr_cols)], id.vars = "name")
ggplot(dot_data, aes(x = variable, y = name, color = value, size = abs(value))) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  theme_minimal() +
  labs(x = "Subtype", y = "Gene", title = "Dot Matrix of Mean Expression")

library(igraph)
motifs <- motifs(ERBB_nodes$g, size = 3) # Find 3-node motifs
barplot(motifs, main = "Motif Frequency")


library(ggplot2)
library(dplyr)
pie_data <- node_metrics %>%
  group_by(community) %>%
  summarise(across(starts_with("mean_expr_"), ~sum(. > 1)))
pie_data_long <- reshape2::melt(pie_data, id.vars = "community")
ggplot(pie_data_long, aes(x = "", y = value, fill = variable)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  facet_wrap(~community) +
  theme_void() +
  labs(title = "Subtype Proportions in Communities")


# Install required packages if needed
# install.packages("BiocManager")
# BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"))

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# Example: genes_of_interest is a vector of gene symbols
 genes_of_interest <- node_metrics$name

# Convert gene symbols to Entrez IDs
gene_df <- bitr(genes_of_interest, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
entrez_genes <- gene_df$ENTREZID

# Perform GO enrichment analysis (Biological Process)
ego <- enrichGO(gene         = entrez_genes,
                OrgDb        = org.Hs.eg.db,
                keyType      = "ENTREZID",
                ont          = "BP",
                pAdjustMethod= "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.2,
                readable     = TRUE)

ego <- pairwise_termsim(ego)
emapplot(ego, showCategory = 30)


CairoPDF("go.pdf", width=15, height=15)
# Visualize as an enrichment map
emapplot(ego, showCategory = 30, label_format = 100, showAll = TRUE)

# Optional: Save the plot
# ggsave("enrichment_map.pdf")
