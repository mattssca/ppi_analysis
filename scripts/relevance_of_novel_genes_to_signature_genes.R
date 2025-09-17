#format data
signature_genes <- uroscanseq_data$signature_genes$ERBB
novel_sig_genes <- subset(erbb_results$sig_genes, !(gene %in% signature_genes)) %>% pull(gene)

##########################################
# Approach 1 - Network Proximity/Shortest Path Analysis
graph_genes <- V(erbb_nodes$g)$name

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

#your vector
rw_scores <- p

#keep only significant genes
rw_scores <- rw_scores[names(rw_scores) %in% novel_sig_genes]

#order by score (descending)
rw_scores <- sort(rw_scores, decreasing = TRUE)

#create data frame for plotting
df_rw <- data.frame(
  Gene = factor(names(rw_scores), levels = names(rw_scores)),
  Score = as.numeric(rw_scores)
)

CairoPDF("analysis_sets/ERBB/figures/big_network/rwr_scores.pdf", width = 20, height = 10)

ggplot(df_rw, aes(x = Gene, y = Score, fill = Score)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_c(option = "viridis") +
  theme_bw() +
  labs(title = "Random Walk with Restart Scores (excluding seeds)",
       x = "Gene", y = "RWR Score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()
##########################################
# Approach 3 - Module/Community Membership

communities <- cluster_louvain(erbb_nodes$g)
membership_vec <- membership(communities)

signature_genes <- c("EGFR", "ERBB2", "ERBB3")
sig_communities <- unique(membership_vec[signature_genes])

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

CairoPDF("analysis_sets/ERBB/figures/big_network/module_membership.pdf", width = 10, height = 20)

#draw plot
ggplot(plot_df, aes(x = module, y = gene, fill = type)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c(Signature = "#0D5EA6", Novel = "#EAA64D")) +
  labs(
    title = "Module Membership of Signature and Novel Signature Genes",
    x = "Module",
    y = "Gene",
    fill = "Gene Type"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

dev.off()


##########################################
# Approach 4 - Functional Similarity (Gene Ontology/Pathways)

#load required packages
library(GOSemSim)
library(org.Hs.eg.db)
library(clusterProfiler)
library(AnnotationDbi)

# Molecular Function (MF): Describes the biochemical activity of a gene product (e.g., kinase activity, DNA binding).
# Cellular Component (CC): Describes where in the cell the gene product is active (e.g., nucleus, plasma membrane).

#map gene symbols to Entrez IDs
gene_map <- bitr(all_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
sig_entrez <- na.omit(gene_map$ENTREZID[gene_map$SYMBOL %in% signature_genes])
novel_entrez <- na.omit(gene_map$ENTREZID[gene_map$SYMBOL %in% novel_sig_genes])
entrez2symbol <- setNames(gene_map$SYMBOL, gene_map$ENTREZID)

#function to calculate overlap for a given ontology
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

#calculate overlaps for BP, MF, CC
overlap_BP <- calc_overlap("BP")
overlap_MF <- calc_overlap("MF")
overlap_CC <- calc_overlap("CC")

#visualization function
plot_overlap <- function(overlap_vec, ontology) {
  df <- data.frame(
    gene = names(overlap_vec),
    overlap = as.numeric(overlap_vec)
  )
  ggplot(head(df, 10), aes(x = reorder(gene, overlap), y = overlap)) +
    geom_col(fill = "#0D5EA6") +
    coord_flip() +
    labs(
      title = paste("GO", ontology, "Term Overlap with Signature Genes"),
      x = "Novel Gene",
      y = paste("Number of Shared", ontology, "Terms")
    ) +
    theme_minimal()
}

#plot for each ontology
plot_overlap(overlap_BP, "BP")
plot_overlap(overlap_MF, "MF")
plot_overlap(overlap_CC, "CC")

#combine into a data frame
df <- data.frame(
  Gene = names(overlap_CC),
  CC = as.numeric(overlap_CC),
  MF = as.numeric(overlap_MF[names(overlap_CC)]),
  BP = as.numeric(overlap_BP[names(overlap_CC)])
)

#if gene order differs, you may want to merge all unique genes:
all_genes <- unique(c(names(overlap_CC), names(overlap_MF), names(overlap_BP)))
df <- data.frame(
  Gene = all_genes,
  CC = as.numeric(overlap_CC[all_genes]),
  MF = as.numeric(overlap_MF[all_genes]),
  BP = as.numeric(overlap_BP[all_genes])
)

#replace NAs with 0 (if a gene is missing in a domain)
df[is.na(df)] <- 0

#convert to long format for ggplot
df_long <- pivot_longer(df, cols = c("CC", "MF", "BP"), names_to = "Domain", values_to = "Overlap")

CairoPDF("analysis_sets/ERBB/figures/big_network/overlap_gene_domain.pdf", width = 20, height = 8)

#plot
ggplot(df_long, aes(x = Gene, y = Overlap, fill = Domain)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  labs(title = "Overlap per Gene and Domain", y = "Overlap", x = "Gene") +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()


########################################
# Combine all!

# Merge all data
df_all <- df_long %>%
  left_join(plot_df, by = c("Gene" = "gene")) %>%
  left_join(df_rw, by = "Gene")

ggplot(df_all, aes(x = Gene, y = Domain)) +
  geom_point(aes(size = Overlap, color = Score, shape = type), alpha = 0.8) +
  scale_color_viridis_c(option = "viridis") +
  scale_size(range = c(3, 20)) +
  theme_bw() +
  labs(title = "Gene Overlap, Domain, and RWR Score",
       color = "RWR Score", size = "Overlap", shape = "Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

library(GGally)
library(dplyr)

# Merge all data
df_all <- df_long %>%
  tidyr::pivot_wider(names_from = Domain, values_from = Overlap) %>%
  left_join(df_rw, by = "Gene") %>%
  left_join(plot_df, by = c("Gene" = "gene"))

# Select and order columns for plotting
df_plot <- df_all %>%
  select(Gene, CC, MF, BP, Score, type)

# Make sure type is a factor
df_plot$type <- as.factor(df_plot$type)

# Parallel coordinates plot
ggparcoord(
  data = df_plot,
  columns = 2:5,           # CC, MF, BP, Score
  groupColumn = 6,         # type
  scale = "uniminmax",     # normalize to [0,1]
  showPoints = TRUE,
  alphaLines = 0.8
) +
  scale_color_manual(values = c(Signature = "#377eb8", Novel = "#e41a1c")) +
  theme_bw() +
  labs(title = "Parallel Coordinates Plot: Overlap & RWR Score by Gene",
       color = "Gene Type")

library(GGally)
library(ggplot2)

# Make sure Gene is a factor
df_plot$Gene <- as.factor(df_plot$Gene)

ggparcoord(
  data = df_plot,
  columns = 2:5,           # CC, MF, BP, Score
  groupColumn = 1,         # Gene
  scale = "uniminmax",
  showPoints = TRUE,
  alphaLines = 0.8
) +
  theme_bw() +
  labs(title = "Parallel Coordinates Plot: Overlap & RWR Score by Gene",
       color = "Gene")
