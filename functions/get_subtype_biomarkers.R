#' @title Identify Subtype-Specific Biomarker Genes for Network Analysis
#'
#' @details This function selects genes that are characteristic of a specific molecular 
#' subtype by identifying genes that are adequately expressed in the target 
#' subtype and show elevated expression compared to other subtypes. The resulting 
#' gene lists are optimized for downstream protein-protein interaction (PPI) 
#' network construction and hub analysis.
#'
#' Algorithm Details:
#' 1. Expression filtering: Genes must have mean expression > min_expr in target subtype
#' 2. Fold change calculation: target_mean / others_mean > fold_change
#' 3. Ranking: Final selection based on highest mean expression in target subtype
#' 4. Missing data: Handled gracefully with na.rm = TRUE in calculations
#' 
#' Parameter Selection Guidelines:
#' Use fold_change = 1 for maximum gene pool (recommended for network analysis)
#' Use fold_change > 1.2 for more stringent subtype-specific selection
#' Adjust min_expr based on your expression data distribution (check quantiles)
#' Set top_n = 200-300 for optimal STRING database mapping success
#'
#' @description
#' The function implements a three-step filtering process:
#' 1. Filter candidate genes to those present in expression data
#' 2. Remove lowly expressed genes (mean expression < min_expr)
#' 3. Calculate fold change vs pooled expression from all other subtypes
#' 4. Keep genes with fold change > threshold (or all if fold_change = 1)
#' 5. If more than top_n genes qualify, select highest expressed ones
#' 
#' This approach maximizes the gene pool for network construction while ensuring 
#' subtype relevance through expression-based selection. Network topology filtering 
#' (hub analysis) will subsequently identify the most functionally important genes.
#'
#' @param expr_target A numeric matrix or data.frame containing gene expression 
#' data for the target subtype. Genes should be in rows, samples in columns.
#' Row names must contain gene symbols.
#' @param expr_others_list A list of numeric matrices or data.frames containing 
#' gene expression data for all other subtypes to compare against. Each element 
#' should have the same structure as expr_target.
#' @param gene_list A character vector of candidate gene symbols to evaluate 
#' (e.g., cell cycle genes, pathway genes). Only genes present in this list 
#' will be considered for selection.
#' @param min_expr Numeric value specifying the minimum mean expression threshold 
#' in the target subtype. Genes with mean expression below this value will be 
#' excluded. Default is 2. Recommended range: 1-3 depending on expression scale.
#' @param fold_change Numeric value specifying the minimum fold change required 
#' vs other subtypes. Set to 1 to include all genes above min_expr (no fold 
#' change filtering). Default is 1.5. Use 1 for maximum gene pool.
#' @param top_n Integer specifying the maximum number of genes to return. If more 
#' genes pass the filtering criteria, the top_n highest expressed genes will be 
#' selected. Default is 150. Recommended range: 100-300.
#'
#' @return A character vector of gene symbols representing potential subtype-specific 
#' biomarkers. Returns an empty character vector if no genes meet the criteria.
#' The genes are ranked by expression level within the target subtype.
#'
get_biomarker_genes <- function(expr_target, 
                                expr_others_list, 
                                gene_list, 
                                min_expr = 2, 
                                fold_change = 1.5, 
                                top_n = 150){
  
  #convert to matrix if it's a data frame (ensure numeric)
  if(is.data.frame(expr_target)) {
    expr_target <- as.matrix(expr_target)
  }
  
  #ensure all others are also matrices
  expr_others_list <- lapply(expr_others_list, function(x) {
    if(is.data.frame(x)) as.matrix(x) else x
  })
  
  #find common genes between gene_list and expression data
  common_genes <- intersect(gene_list, rownames(expr_target))
  
  if(length(common_genes) == 0) {
    warning("No genes from gene_list found in expression data")
    return(character(0))
  }
  
  #step 1: Filter by minimum expression
  mean_expr_target <- rowMeans(expr_target[common_genes, , drop = FALSE], na.rm = TRUE)
  expressed_genes <- names(mean_expr_target[mean_expr_target > min_expr & !is.na(mean_expr_target)])
  
  if(length(expressed_genes) == 0) {
    warning("No genes meet minimum expression threshold")
    return(character(0))
  }
  
  #step 2: Get fold change vs other subtypes
  fc_genes <- c()
  for(gene in expressed_genes) {
    target_mean <- mean(as.numeric(expr_target[gene, ]), na.rm = TRUE)
    
    #calculate mean across all other subtypes
    other_values <- c()
    for(other_expr in expr_others_list) {
      if(gene %in% rownames(other_expr)) {
        other_values <- c(other_values, as.numeric(other_expr[gene, ]))
      }
    }
    
    if(length(other_values) > 0) {
      others_mean <- mean(other_values, na.rm = TRUE)
      
      #check for valid values before comparison
      if(!is.na(target_mean) && !is.na(others_mean) && others_mean > 0 && target_mean/others_mean > fold_change) {
        fc_genes <- c(fc_genes, gene)
      }
    }
  }
  
  #step 3: Rank by expression and take top N
  if(length(fc_genes) > top_n) {
    fc_expr <- mean_expr_target[fc_genes]
    fc_genes <- names(sort(fc_expr, decreasing = TRUE)[1:top_n])
  }
  
  return(fc_genes)
}
