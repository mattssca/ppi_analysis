#load libraries
library(DOSE)
library(clusterProfiler)
library(tidyr)

#convert gene symbols to Entrez IDs
library(org.Hs.eg.db)
entrez_ids <- bitr(novel_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID

#disease Ontology enrichment
do_enrich <- enrichDO(gene = entrez_ids, ont = "DO", readable = TRUE)

#view results for bladder cancer
subset(do_enrich@result, grepl("bladder", Description, ignore.case = TRUE))

#select relevant columns and split geneID into a vector of gene symbols
output <- do_enrich[, c("ID", "Description", "geneID", "Count", "p.adjust")]

#subset data to bladder related genes
bladder_hits = output %>% 
  filter(Description %in% c("urinary bladder cancer", "bladder carcinoma", "bladder disease"))

#split geneID into a list of gene symbols for each row
bladder_hits$GeneSymbols <- strsplit(as.character(bladder_hits$geneID), "/")

#print the genes for each description
for (i in seq_len(nrow(bladder_hits))) {
  cat("Description:", bladder_hits$Description[i], "\n")
  cat("Genes:", paste(bladder_hits$GeneSymbols[[i]], collapse = ", "), "\n\n")
}

#split geneID into a list of gene symbols
bladder_hits$Gene <- strsplit(as.character(bladder_hits$geneID), "/")

#unnest so each row is one gene
bladder_long <- bladder_hits %>%
  select(ID, Description, Count, p.adjust, Gene) %>%
  unnest(Gene)

#view the result
head(bladder_long)
head(novel_sig_genes)
head(signature_genes)

#genes annotated as bladder-related
bladder_genes <- unique(bladder_long$Gene)

#annotate each novel gene as "Bladder-related" or "Not annotated"
novel_annot_df_condensed <- data.frame(
  gene = novel_sig_genes,
  description_condensed = ifelse(novel_sig_genes %in% bladder_genes, "Bladder-related", "Not annotated")
)

head(novel_annot_df_condensed)

#for each gene, get all associated bladder-related descriptions (or "Not annotated")
novel_annot_df <- data.frame(
  gene = novel_sig_genes,
  description = sapply(novel_sig_genes, function(g) {
    descs <- unique(bladder_long$Description[bladder_long$Gene == g])
    if (length(descs) == 0) {
      "Not annotated"
    } else {
      paste(descs, collapse = "; ")
    }
  }),
  stringsAsFactors = FALSE
)

#set levels
novel_annot_df$description <- factor(novel_annot_df$description)
head(novel_annot_df)

#join them
bladder_genes_association = novel_annot_df %>% 
  left_join(novel_annot_df_condensed, by = "gene")


#make sure description_condensed is a factor with desired order
bladder_genes_association$description_condensed <- factor(
  bladder_genes_association$description_condensed,
  levels = c("Bladder-related", "Not annotated")
)

#barplot: count of genes by annotation status
ggplot(bladder_genes_association, aes(x = description_condensed, fill = description_condensed)) +
  geom_bar() +
  scale_fill_manual(values = c("Bladder-related" = "#FAA533", "Not annotated" = "#44444E")) +
  labs(
    title = "Novel Significant Genes: Bladder Cancer Annotation",
    x = "Annotation Status",
    y = "Gene Count",
    fill = "Annotation"
  ) +
  theme_minimal()


# Merge hub_score into plot_df
plot_df <- merge(
  plot_df,
  erbb_nodes$node_metrics[, c("name", "hub_score")],
  by.x = "gene",
  by.y = "name",
  all.x = TRUE
)

# For plotting, calculate -log10(padj)
plot_df$neglog10_padj <- -log10(plot_df$padj)

# Order genes by effect size or significance for plotting
plot_df$gene <- factor(plot_df$gene, levels = plot_df$gene[order(plot_df$neglog10_padj, decreasing = TRUE)])

# Rank genes by significance or effect size
plot_df <- plot_df[order(plot_df$neglog10_padj, decreasing = TRUE), ]
plot_df$rank <- seq_len(nrow(plot_df))

# Add "Seed" as a level if not present
if (!"Seed" %in% levels(plot_df$description)) {
  plot_df$description <- factor(plot_df$description, levels = c(levels(plot_df$description), "Seed"))
}
if (!"Seed" %in% levels(plot_df$description_condensed)) {
  plot_df$description_condensed <- factor(plot_df$description_condensed, levels = c(levels(plot_df$description_condensed), "Seed"))
}

# Now assign
plot_df$description[plot_df$gene %in% signature_genes] <- "Seed"
plot_df$description_condensed[plot_df$gene %in% signature_genes] <- "Seed"

ggplot(plot_df, aes(x = rank, y = neglog10_padj, color = description_condensed, label = gene, size = hub_score)) +
  geom_point(alpha = 0.9) +
  geom_text_repel(size = 3, max.overlaps = 20) +
  scale_color_manual(values = c("Bladder-related" = "#FAA533", "Not annotated" = "#44444E", "Seed" = "#77BEF0")) +
  scale_size_continuous(range = c(0.5, 5)) +  # adjust dot size range as needed
  labs(
    x = "Gene Rank (by significance)",
    y = expression(-log[10]("adj. p-value")),
    color = "Bladder Annotation",
    size = "Hub Score",
    title = "Rainfall Plot of Gene Significance"
  ) +
  theme_minimal()
library(ggplot2)

# Order genes by significance
plot_df <- plot_df[order(plot_df$neglog10_padj, decreasing = TRUE), ]
plot_df$gene <- factor(plot_df$gene, levels = plot_df$gene)

ggplot(plot_df, aes(x = neglog10_padj, y = gene, color = description_condensed, size = hub_score)) +
  geom_segment(aes(x = 0, xend = neglog10_padj, y = gene, yend = gene), size = 0.2, color = "black") +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("Bladder-related" = "#FAA533", "Not annotated" = "#44444E", "Seed" = "#77BEF0")) +
  scale_size_continuous(range = c(0.1, 6)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
  labs(
    x = expression(-log[10]("adj. p-value")),
    y = "Gene (ranked by significance)",
    color = "Bladder Annotation",
    size = "Hub Score",
    title = "Lollipop Plot of Gene Significance"
  ) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 6, face = "bold"))
