# Test if hub scores differ significantly between subtypes
library(broom)
library(networkD3)
library(igraph)
library(clusterProfiler)
library(org.Hs.eg.db)
library(UpSetR)


# Test if subtype-specific hubs are enriched for specific pathways
for(subtype in names(subtype_hubs)) {
  top_genes <- subtype_hubs[[subtype]]$gene[1:10]
  
  # GO enrichment
  enrichment <- enrichGO(gene = top_genes,
                         OrgDb = org.Hs.eg.db,
                         keyType = "SYMBOL",
                         ont = "BP",
                         pAdjustMethod = "BH")
  
  print(paste("Enrichment for", subtype))
  print(head(enrichment@result))
}

# For a single subtype
enrichment_df <- enrichment@result %>%
  slice_head(n = 10) %>%
  mutate(Description = reorder(Description, -log10(p.adjust)))

ggplot(enrichment_df, aes(x = -log10(p.adjust), y = Description)) +
  geom_col(fill = "steelblue", alpha = 0.7) +
  labs(title = "GO Enrichment Analysis",
       x = "-log10(Adjusted P-value)",
       y = "GO Terms") +
  theme_minimal()

# Initialize empty dataframe
library(stringr)
library(dplyr)

# Define the Lund colors
lund_colors <- c(
  "uro" = "#3cb44b",
  "gu" = "#4363d8", 
  "basq" = "#CD2626",
  "mes" = "#f58231",
  "scne" = "#A020F0"
)

all_enrichments <- data.frame()

for(subtype in names(subtype_hubs)) {
  top_genes <- subtype_hubs[[subtype]]$gene[1:10]
  
  # Complete enrichGO function call with all required parameters
  enrich <- enrichGO(gene = top_genes,
                     OrgDb = org.Hs.eg.db,
                     keyType = "SYMBOL",
                     ont = "BP",  # Biological Process
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.2)
  
  if(nrow(enrich@result) > 0) {
    temp_df <- enrich@result %>%
      slice_head(n = 5) %>%
      mutate(subtype = subtype,
             # Simplify common GO term patterns
             Description_simple = Description %>%
               str_replace_all("regulation of ", "reg. ") %>%
               str_replace_all("positive regulation of ", "+reg. ") %>%
               str_replace_all("negative regulation of ", "-reg. ") %>%
               str_replace_all("biological process", "process") %>%
               str_replace_all("cellular component", "component") %>%
               str_replace_all("molecular function", "function"),
             # Wrap text to multiple lines
             Description_wrapped = str_wrap(Description_simple, width = 30))
    all_enrichments <- bind_rows(all_enrichments, temp_df)
  }
}  # <- Missing closing brace was here

ggplot(all_enrichments, aes(x = -log10(p.adjust), y = Description_simple, fill = subtype)) +
  geom_col(position = "dodge") +
  facet_wrap(~subtype, scales = "free_y") +
  scale_fill_manual(values = lund_colors) +  # <- Added Lund colors
  labs(title = "GO Enrichment by Subtype",
       y = "GO Terms") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

# Show which genes overlap between enriched terms
gene_sets <- strsplit(enrichment@result$geneID[1:10], "/")
names(gene_sets) <- enrichment@result$Description[1:10]

upset(fromList(gene_sets), 
      order.by = "freq",
      nsets = 10,
      main.bar.color = "steelblue",
      sets.bar.color = "darkred")

# Heatmap of enrichment across subtypes
library(pheatmap)

enrichment_matrix <- all_enrichments %>%
  dplyr::select(Description, subtype, p.adjust) %>%
  dplyr::group_by(Description) %>%
  dplyr::mutate(
    unique_desc = ifelse(n() > 1, 
                         paste0(Description, " (", subtype, ")"), 
                         Description)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(neg_log_p = -log10(p.adjust)) %>%
  dplyr::select(unique_desc, subtype, neg_log_p) %>%
  tidyr::pivot_wider(names_from = subtype, 
                     values_from = neg_log_p, 
                     values_fill = 0) %>%
  tibble::column_to_rownames("unique_desc") %>%
  as.matrix()

pheatmap(enrichment_matrix,
         color = colorRampPalette(c("white", "red"))(100),
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         main = "GO Enrichment Heatmap\n(-log10 adjusted p-value)",
         fontsize_row = 8)
