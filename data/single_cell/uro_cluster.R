# Complete Urothelial Binary Classification Analysis
library(Seurat)
library(dplyr)

# Function to perform binary urothelial classification
classify_urothelial_genes <- function(seurat_object, 
                                      gene_list, 
                                      urothelial_cluster = 0,
                                      detection_threshold = 0.05,
                                      expression_threshold = 0.1,
                                      strict_detection = 0.1,
                                      strict_expression = 0.2,
                                      lenient_detection = 0.01,
                                      lenient_expression = 0.05) {
  
  # Ensure the object has the required cluster information
  if(!"seurat_clusters" %in% colnames(seurat_object@meta.data)) {
    stop("Seurat object must have 'seurat_clusters' in meta.data. Run FindClusters() first.")
  }
  
  cat("Performing binary urothelial classification...\n")
  cat("Urothelial cluster:", urothelial_cluster, "\n")
  cat("Total genes to analyze:", length(gene_list), "\n")
  
  # Initialize results dataframe
  urothelial_detection <- data.frame()
  
  # Loop through genes
  for(gene in gene_list) {
    if(gene %in% rownames(seurat_object)) {
      # Get expression in urothelial cells only
      gene_expr <- GetAssayData(seurat_object, slot = "data")[gene, ]
      cluster_ids <- seurat_object@meta.data$seurat_clusters
      
      urothelial_expr <- gene_expr[cluster_ids == urothelial_cluster]
      
      detection_rate <- sum(urothelial_expr > 0) / length(urothelial_expr)
      mean_expr <- mean(urothelial_expr)
      max_expr <- max(urothelial_expr)
      
      urothelial_detection <- rbind(urothelial_detection, data.frame(
        Gene = gene,
        Mean_Expression = mean_expr,
        Detection_Rate = detection_rate,
        Max_Expression = max_expr,
        Cells_Expressing = sum(urothelial_expr > 0),
        Total_Urothelial_Cells = length(urothelial_expr),
        stringsAsFactors = FALSE
      ))
    } else {
      # Gene not found in dataset
      urothelial_detection <- rbind(urothelial_detection, data.frame(
        Gene = gene,
        Mean_Expression = 0,
        Detection_Rate = 0,
        Max_Expression = 0,
        Cells_Expressing = 0,
        Total_Urothelial_Cells = sum(seurat_object@meta.data$seurat_clusters == urothelial_cluster),
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Standard binary classification
  urothelial_detection$Classification <- ifelse(
    urothelial_detection$Detection_Rate >= detection_threshold & 
      urothelial_detection$Mean_Expression >= expression_threshold,
    "Urothelial",
    "Non-urothelial"
  )
  
  # Strict criteria
  urothelial_detection$Strict_Classification <- ifelse(
    urothelial_detection$Detection_Rate >= strict_detection & 
      urothelial_detection$Mean_Expression >= strict_expression,
    "Urothelial",
    "Non-urothelial"
  )
  
  # Lenient criteria  
  urothelial_detection$Lenient_Classification <- ifelse(
    urothelial_detection$Detection_Rate >= lenient_detection | 
      urothelial_detection$Mean_Expression >= lenient_expression,
    "Urothelial",
    "Non-urothelial"
  )
  
  # Sort by mean expression
  urothelial_detection <- urothelial_detection[order(urothelial_detection$Mean_Expression, decreasing = TRUE), ]
  
  # Print summary
  cat("\n=== CLASSIFICATION RESULTS ===\n")
  cat("Standard criteria (detection >=", detection_threshold, "& expression >=", expression_threshold, "):\n")
  print(table(urothelial_detection$Classification))
  
  cat("\nStrict criteria (detection >=", strict_detection, "& expression >=", strict_expression, "):\n")
  print(table(urothelial_detection$Strict_Classification))
  
  cat("\nLenient criteria (detection >=", lenient_detection, "OR expression >=", lenient_expression, "):\n")
  print(table(urothelial_detection$Lenient_Classification))
  
  # Show urothelial genes
  urothelial_genes <- urothelial_detection[urothelial_detection$Classification == "Urothelial", ]
  cat("\nGenes classified as UROTHELIAL (standard criteria):\n")
  print(urothelial_genes[, c("Gene", "Mean_Expression", "Detection_Rate", "Cells_Expressing")])
  
  cat("\nSummary:\n")
  cat("Total genes analyzed:", nrow(urothelial_detection), "\n")
  cat("Genes found in dataset:", sum(urothelial_detection$Mean_Expression > 0 | urothelial_detection$Detection_Rate > 0), "\n")
  cat("Standard urothelial genes:", nrow(urothelial_genes), "\n")
  cat("Percentage urothelial (standard):", round(nrow(urothelial_genes)/nrow(urothelial_detection)*100, 1), "%\n")
  
  return(urothelial_detection)
}

# ===============================================================
# MAIN ANALYSIS WORKFLOW
# ===============================================================

# Step 1: Load and prepare data
bladder_combined <- bladder_normal  # or bladder_cancer

# Step 2: Standard single-cell processing
cat("Starting single-cell processing...\n")
bladder_combined <- NormalizeData(bladder_combined)
bladder_combined <- FindVariableFeatures(bladder_combined)
bladder_combined <- ScaleData(bladder_combined)
bladder_combined <- RunPCA(bladder_combined)

# Step 3: Clustering
cat("Performing clustering...\n")
n_pcs <- min(15, ncol(Embeddings(bladder_combined, "pca")))
bladder_combined <- FindNeighbors(bladder_combined, dims = 1:n_pcs)
bladder_combined <- FindClusters(bladder_combined, resolution = 0.5)
bladder_combined <- RunUMAP(bladder_combined, dims = 1:n_pcs)
bladder_combined <- JoinLayers(bladder_combined)

# Step 4: Visualize clusters
cat("Visualizing clusters...\n")
DimPlot(bladder_combined, reduction = "umap", label = TRUE, label.size = 6)

# Step 5: Validate urothelial cluster identity
cat("Validating urothelial cluster identity...\n")
urothelial_markers <- c("KRT20", "UPK1A", "UPK2", "UPK3A", "GATA3", "FOXA1", "KRT7", "KRT8")
available_uro_markers <- intersect(urothelial_markers, rownames(bladder_combined))

if(length(available_uro_markers) > 0) {
  cat("Available urothelial markers:", paste(available_uro_markers, collapse = ", "), "\n")
  
  # Plot markers
  FeaturePlot(bladder_combined, features = available_uro_markers[1:min(4, length(available_uro_markers))], ncol = 2)
  
  # Check expression by cluster
  marker_expr <- AverageExpression(bladder_combined, 
                                   features = available_uro_markers,
                                   group.by = "seurat_clusters")$RNA
  cat("Urothelial marker expression by cluster:\n")
  print(marker_expr)
  
  # Identify cluster with highest urothelial marker expression
  cluster_scores <- colMeans(marker_expr, na.rm = TRUE)
  urothelial_cluster_id <- as.numeric(names(which.max(cluster_scores)))
  cat("Predicted urothelial cluster:", urothelial_cluster_id, "\n")
} else {
  cat("No urothelial markers found in dataset, assuming cluster 0\n")
  urothelial_cluster_id <- 0
}

# Step 6: Perform binary classification of your ERBB genes
cat("\nStarting ERBB gene classification...\n")
results <- classify_urothelial_genes(
  seurat_object = bladder_combined,
  gene_list = erbb_genes,
  urothelial_cluster = urothelial_cluster_id,
  detection_threshold = 0.05,
  expression_threshold = 0.1
)

# Step 7: Additional analysis options
cat("\n=== ADDITIONAL ANALYSIS OPTIONS ===\n")

# Option A: Strict criteria for high-confidence genes
cat("Running strict criteria analysis...\n")
strict_results <- classify_urothelial_genes(
  seurat_object = bladder_combined,
  gene_list = erbb_genes,
  urothelial_cluster = urothelial_cluster_id,
  detection_threshold = 0.1,
  expression_threshold = 0.2
)

# Option B: Merge with your protein atlas data if available
if(exists("erbb_extended")) {
  cat("Merging with protein atlas data...\n")
  merged_results <- merge(results, erbb_extended, by.x = "Gene", by.y = "gene_symbol", all.x = TRUE)
  
  # Focus on protein atlas urothelium-active genes
  uro_active_genes <- merged_results[merged_results$protein_atlas == "uro_cells" & !is.na(merged_results$protein_atlas), ]
  if(nrow(uro_active_genes) > 0) {
    cat("Genes that are both scRNA-seq urothelial AND protein atlas urothelium-active:\n")
    uro_validated <- uro_active_genes[uro_active_genes$Classification == "Urothelial", ]
    print(uro_validated[, c("Gene", "Mean_Expression", "Detection_Rate", "Classification", "protein_atlas")])
  }
}

# Option C: Visualize top urothelial genes
top_urothelial_genes <- head(results[results$Classification == "Urothelial", ]$Gene, 6)
if(length(top_urothelial_genes) > 0) {
  cat("Visualizing top urothelial ERBB genes...\n")
  FeaturePlot(bladder_combined, features = top_urothelial_genes, ncol = 3)
  VlnPlot(bladder_combined, features = top_urothelial_genes[1:4], ncol = 2, group.by = "seurat_clusters")
}

# Step 8: Save results
cat("Analysis complete! Results saved to 'results' object\n")
cat("Access results with: View(results)\n")
cat("Access merged results (if available) with: View(merged_results)\n")

# Print final summary
cat("\n=== FINAL SUMMARY ===\n")
cat("Total ERBB genes analyzed:", nrow(results), "\n")
cat("Standard urothelial genes:", sum(results$Classification == "Urothelial"), "\n")
cat("Strict urothelial genes:", sum(results$Strict_Classification == "Urothelial"), "\n")
cat("Lenient urothelial genes:", sum(results$Lenient_Classification == "Urothelial"), "\n")