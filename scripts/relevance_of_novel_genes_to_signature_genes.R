##########################################
# Approach 1 - Network Proximity/Shortest Path Analysis
graph_genes <- V(erbb_nodes$g)$name
signature_genes <- c("EGFR", "ERBB2", "ERBB3")
novel_sig_genes <- subset(sig_genes, !(gene %in% signature_genes)) %>% pull(gene)

novel_sig_genes_in_graph <- intersect(novel_sig_genes, graph_genes)
signature_genes_in_graph <- intersect(signature_genes, graph_genes)

distances_to_signature <- sapply(novel_sig_genes, function(gene) {
  min(distances(erbb_nodes$g, v = gene, to = signature_genes))
})


##########################################
# Appraoch 2 - Network Propagation/Random Walk with Restart

#make sure vertex names are gene symbols
adjmat <- igraph::as_adjacency_matrix(erbb_nodes$g, sparse = FALSE)

#all genes in the network
all_genes <- rownames(adjmat)

#normalize adjacency matrix by columns
normalize <- function(mat) {
  sweep(mat, 2, colSums(mat), "/")
}
A <- normalize(adjmat)

#seed vector: 1 for signature genes, 0 for others
all_genes <- rownames(A)
p0 <- as.numeric(all_genes %in% signature_genes)
names(p0) <- all_genes

#parameters
restart <- 0.7
tol <- 1e-10
p <- p0
repeat {
  p_new <- (1 - restart) * (A %*% p) + restart * p0
  if (sum(abs(p_new - p)) < tol) break
  p <- p_new
}

#if p is a matrix, convert to vector
if (is.matrix(p)) p <- as.vector(p)
names(p) <- rownames(adjmat)

#sort and select top 20
top_genes <- sort(p, decreasing = TRUE)[1:20]

barplot(top_genes, las = 2, col = "steelblue",
        main = "Top 20 Network Propagation Scores",
        ylab = "Propagation Score")

##########################################
# Approach 3 - Module/Community Membership

communities <- cluster_louvain(erbb_nodes$g)
membership_vec <- membership(communities)

signature_genes <- c("EGFR", "ERBB2", "ERBB3")
sig_communities <- unique(membership_vec[signature_genes])

#get your novel genes (excluding signature genes)
novel_sig_genes <- setdiff(sig_genes$gene, signature_genes)

#which novel genes are in the same community as any signature gene?
novel_in_sig_module <- novel_sig_genes[membership_vec[novel_sig_genes] %in% sig_communities]

#optionally, create a data frame with module info
module_df <- data.frame(
  gene = novel_sig_genes,
  module = membership_vec[novel_sig_genes],
  in_signature_module = membership_vec[novel_sig_genes] %in% sig_communities
)

#combine signature and novel genes
all_genes <- c(signature_genes, novel_sig_genes)
gene_type <- c(rep("Signature", length(signature_genes)), rep("Novel", length(novel_sig_genes)))
gene_module <- membership_vec[all_genes]

#prepare plot data
plot_df <- data.frame(
  gene = all_genes,
  module = as.factor(gene_module),
  type = gene_type
)

#perserve order
plot_df$gene <- factor(plot_df$gene, levels = unique(plot_df$gene))

#draw plot
ggplot(plot_df, aes(x = module, y = gene, fill = type)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c(Signature = "red", Novel = "skyblue")) +
  labs(
    title = "Module Membership of Signature and Novel Signature Genes",
    x = "Module",
    y = "Gene",
    fill = "Gene Type"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))


##########################################
# Approach 4 - Functional Similarity (Gene Ontology/Pathways)

# Load required packages
library(GOSemSim)
library(org.Hs.eg.db)
library(clusterProfiler)
library(AnnotationDbi)
library(ggplot2)

# Molecular Function (MF): Describes the biochemical activity of a gene product (e.g., kinase activity, DNA binding).
# Cellular Component (CC): Describes where in the cell the gene product is active (e.g., nucleus, plasma membrane).

# Example gene sets
signature_genes <- c("EGFR", "ERBB2", "ERBB3")
novel_sig_genes <- c("SRC", "STAT3", "PDGFRB", "TP53", "JAK2", "LYN")

# Map gene symbols to Entrez IDs
all_genes <- unique(c(signature_genes, novel_sig_genes))
gene_map <- bitr(all_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
sig_entrez <- na.omit(gene_map$ENTREZID[gene_map$SYMBOL %in% signature_genes])
novel_entrez <- na.omit(gene_map$ENTREZID[gene_map$SYMBOL %in% novel_sig_genes])
entrez2symbol <- setNames(gene_map$SYMBOL, gene_map$ENTREZID)

# Function to calculate overlap for a given ontology
calc_overlap <- function(ontology) {
  sig_go_terms <- AnnotationDbi::select(org.Hs.eg.db, keys = sig_entrez, columns = "GO", keytype = "ENTREZID")
  sig_go_terms <- subset(sig_go_terms, ONTOLOGY == ontology)
  novel_go_terms <- AnnotationDbi::select(org.Hs.eg.db, keys = novel_entrez, columns = "GO", keytype = "ENTREZID")
  novel_go_terms <- subset(novel_go_terms, ONTOLOGY == ontology)
  sig_terms <- unique(sig_go_terms$GO)
  overlap <- sapply(unique(novel_go_terms$ENTREZID), function(novel_id) {
    novel_terms <- novel_go_terms$GO[novel_go_terms$ENTREZID == novel_id]
    sum(novel_terms %in% sig_terms)
  })
  names(overlap) <- entrez2symbol[names(overlap)]
  sort(overlap, decreasing = TRUE)
}

# Calculate overlaps for BP, MF, CC
overlap_BP <- calc_overlap("BP")
overlap_MF <- calc_overlap("MF")
overlap_CC <- calc_overlap("CC")

# Visualization function
plot_overlap <- function(overlap_vec, ontology) {
  df <- data.frame(
    gene = names(overlap_vec),
    overlap = as.numeric(overlap_vec)
  )
  ggplot(head(df, 10), aes(x = reorder(gene, overlap), y = overlap)) +
    geom_col(fill = "salmon") +
    coord_flip() +
    labs(
      title = paste("GO", ontology, "Term Overlap with Signature Genes"),
      x = "Novel Gene",
      y = paste("Number of Shared", ontology, "Terms")
    ) +
    theme_minimal()
}

# Plot for each ontology
plot_overlap(overlap_BP, "BP")
plot_overlap(overlap_MF, "MF")
plot_overlap(overlap_CC, "CC")

