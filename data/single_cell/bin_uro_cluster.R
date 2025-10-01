library(Seurat)
library(dplyr)

seurat_object = bladder_cancer_processed

#define known urothelial markers
urothelial_markers <- c("GATA3", "FOXA1", "KRT7", "KRT8", "PPARG", "FOXA1", 
                        "ELF3", "FGFR3", "CCND1", "TP63", "KRT5", "KRT14", 
                        "EPCAM", "CDH1")

#fix data layers (prevents errors)
seurat_object <- JoinLayers(seurat_object)

#step 2: check which urothelial markers are available in the data
available_uro_markers <- intersect(urothelial_markers, rownames(seurat_object))
cat("Available urothelial markers:", paste(available_uro_markers, collapse = ", "), "\n")

#step 3: find which cluster expresses urothelial markers most
if(length(available_uro_markers) > 0) {
  marker_expr <- AverageExpression(seurat_object, 
                                   features = available_uro_markers,
                                   group.by = "seurat_clusters")$RNA
  
  cat("Urothelial marker expression by cluster:\n")
  print(marker_expr)
  
  #check what the cluster names look like
  cat("Cluster names:", colnames(marker_expr), "\n")
  
  #calculate cluster scores and find max
  cluster_scores <- colMeans(marker_expr, na.rm = TRUE)
  cat("Cluster scores:", cluster_scores, "\n")
  
  #get the cluster with highest score (keep as character)
  max_cluster_name <- names(which.max(cluster_scores))
  cat("Max cluster name:", max_cluster_name, "\n")
  
  #keep as character (don't convert to numeric since they're "g0", "g1")
  urothelial_cluster <- max_cluster_name
  
  cat("Urothelial cluster identified as:", urothelial_cluster, "\n")
} else {
  # If no markers found, assume first cluster
  unique_clusters <- unique(seurat_object@meta.data$seurat_clusters)
  urothelial_cluster <- unique_clusters[1]
  cat("No urothelial markers found, assuming cluster:", urothelial_cluster, "\n")
}

# Verify the cluster assignment
cat("Final urothelial cluster:", urothelial_cluster, "\n")
cat("This cluster has higher expression of urothelial markers\n")

#step 4: Get cluster assignments for all cells
cluster_ids <- seurat_object@meta.data$seurat_clusters

#step 5: Set up classification parameters
detection_threshold <- 0.05  # 5% of cells must express the gene
expression_threshold <- 0.1   # Minimum average expression level

#fix the cluster assignment
#g0 had the highest urothelial marker expression, use cluster 0
urothelial_cluster_id <- 0
cat("Using cluster ID:", urothelial_cluster_id, "for urothelial cells\n")

#verify this works
cat("Number of cells in cluster 0:", sum(cluster_ids == 0), "\n")
cat("Number of cells in cluster 1:", sum(cluster_ids == 1), "\n")

#step 6: Binary classify each ERBB gene
classification_results <- data.frame()

for(gene in erbb_genes) {
  if(gene %in% rownames(seurat_object)) {
    
    #get expression for this gene
    gene_expr <- GetAssayData(seurat_object, slot = "data")[gene, ]
    
    #get expression only in urothelial cells (cluster 0)
    urothelial_expr <- gene_expr[cluster_ids == urothelial_cluster_id]
    
    #calculate statistics
    mean_expr <- mean(urothelial_expr)
    detection_rate <- sum(urothelial_expr > 0) / length(urothelial_expr)
    max_expr <- max(urothelial_expr)
    cells_expressing <- sum(urothelial_expr > 0)
    
    #binary classification
    classification <- ifelse(
      detection_rate >= detection_threshold & mean_expr >= expression_threshold,
      "Urothelial",
      "Non-urothelial"
    )
    
    #store results
    classification_results <- rbind(classification_results, data.frame(
      Gene = gene,
      Mean_Expression = mean_expr,
      Detection_Rate = detection_rate,
      Max_Expression = max_expr,
      Cells_Expressing = cells_expressing,
      Total_Urothelial_Cells = length(urothelial_expr),
      Classification = classification,
      stringsAsFactors = FALSE
    ))
    
  } else {
    # Gene not found in dataset
    classification_results <- rbind(classification_results, data.frame(
      Gene = gene,
      Mean_Expression = 0,
      Detection_Rate = 0,
      Max_Expression = 0,
      Cells_Expressing = 0,
      Total_Urothelial_Cells = sum(cluster_ids == urothelial_cluster_id),
      Classification = "Non-urothelial",
      stringsAsFactors = FALSE
    ))
  }
}

#step 7: sort results by expression level
classification_results <- classification_results[order(classification_results$Mean_Expression, decreasing = TRUE), ]

#step 8: view results
cat("\n=== CLASSIFICATION SUMMARY ===\n")
table(classification_results$Classification)

cat("\nGenes classified as UROTHELIAL:\n")
urothelial_genes <- classification_results[classification_results$Classification == "Urothelial", ]
print(urothelial_genes[, c("Gene", "Mean_Expression", "Detection_Rate", "Cells_Expressing")])

cat("\nTotal ERBB genes analyzed:", nrow(classification_results), "\n")
cat("Urothelial genes:", nrow(urothelial_genes), "\n")
cat("Percentage urothelial:", round(nrow(urothelial_genes)/nrow(classification_results)*100, 1), "%\n")

# Step 9: View full results
View(classification_results)
