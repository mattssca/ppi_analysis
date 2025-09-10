#' Run ANOVA and Post-hoc Analysis for PPI Genes by Subtype
#'
#' Performs one-way ANOVA for each gene in a set of extended PPI genes to test for differences in expression across molecular subtypes. 
#' Significant genes are identified based on adjusted p-value and effect size. Optionally, Tukey's HSD post-hoc tests are performed for significant genes.
#' Generates boxplots of expression by subtype for all significant genes.
#'
#' @param expr_data A numeric matrix or data frame of gene expression values (genes x samples). Row names are gene symbols, column names are sample IDs.
#' @param subtype_vector A named vector or factor of subtype assignments for each sample. Names must match column names of \code{expr_data}.
#' @param extended_genes A character vector of gene symbols to test (must be present in \code{expr_data}).
#' @param subtype_class Character. Either \code{"5_class"} or \code{"7_class"}; determines the order of subtypes for plotting.
#' @param sig_p Numeric. Adjusted p-value threshold for significance (default: 0.01).
#' @param min_diff Numeric. Minimum difference in mean expression between subtypes required for significance (default: 1).
#'
#' @return A list with the following elements:
#'   \item{sig_genes}{A data frame of significant genes with p-values, adjusted p-values, and effect sizes.}
#'   \item{anova_results}{A data frame of ANOVA results for all tested genes.}
#'   \item{posthoc_results}{A named list of Tukey HSD post-hoc results for each significant gene.}
#'   \item{plot}{A patchwork object containing boxplots of expression by subtype for all significant genes.}
#'
#' @details
#' For each gene, a one-way ANOVA is performed to test for differences in mean expression across subtypes. 
#' P-values are adjusted for multiple testing using the Bonferroni method. Genes passing the significance and effect size thresholds are subjected to Tukey's HSD post-hoc test.
#' Boxplots are generated for all significant genes, colored by subtype.
#'
#' @import ggplot2
#' @import patchwork
#' @export
#'
#' @examples
#' \dontrun{
#' results <- run_anova_ppi(expr_data, subtype_vector, extended_genes, subtype_class = "5_class")
#' }
run_anova_ppi = function(expr_data, 
                         subtype_vector, 
                         extended_genes,
                         node_metrics, 
                         subtype_class = "5_class", 
                         sig_p = 0.01, 
                         min_diff = 1){
  
  #load packages
  library(ggplot2)
  library(patchwork)
  
  #subset expression to genes of interest
  expr_sub <- expr_data[extended_genes, , drop = FALSE]
  
  #run ANOVA
  anova_results <- lapply(rownames(expr_sub), function(extended_genes) {
    expr_vec <- as.numeric(expr_sub[extended_genes, ])
    df <- data.frame(expr = expr_vec, subtype = factor(subtype_vector[colnames(expr_sub)]))
    fit <- aov(expr ~ subtype, data = df)
    pval <- summary(fit)[[1]][["Pr(>F)"]][1]
    means <- tapply(df$expr, df$subtype, mean)
    max_diff <- max(means) - min(means)
    data.frame(gene = extended_genes, pval = pval, max_diff = max_diff)
  })
  
  #add p value correction method
  anova_df <- do.call(rbind, anova_results)
  anova_df$padj <- p.adjust(anova_df$pval, method = "bonferroni")
  
  #filter for both significance and effect size
  sig_genes <- subset(anova_df, padj < sig_p & max_diff > min_diff)
  
  #run post hoc
  posthoc_results <- lapply(sig_genes$gene, function(extended_genes) {
    expr_vec <- as.numeric(expr_sub[extended_genes, ])
    df <- data.frame(expr = expr_vec, subtype = factor(subtype_vector[colnames(expr_sub)]))
    fit <- aov(expr ~ subtype, data = df)
    tukey <- TukeyHSD(fit)
    sig_pairs <- subset(as.data.frame(tukey$subtype), `p adj` < 0.05)
    list(gene = extended_genes, tukey = tukey, sig_pairs = sig_pairs)
    })
  names(posthoc_results) <- sig_genes$gene

  #plots
  if(subtype_class == "5_class"){
    subtype_order <- c("Uro", "GU", "BaSq")
  }else if(subtype_class == "7_class"){
    subtype_order <- c("UroA", "UroB", "UroC", "GU", "BaSq")
  }else{
    stop("Please provide a valid subtype class (5_class or 7_class)...")
  }

  plots <- list()
  
  for (i in seq_along(sig_genes$gene)) {
    gene <- sig_genes$gene[i]
    expr_vec <- as.numeric(expr_sub[gene, ])
    df <- data.frame(expr = expr_vec, subtype = factor(subtype_vector[colnames(expr_sub)], levels = subtype_order))    
    padj <- signif(sig_genes$padj[i], 2)
    p <- ggplot(df, aes(x = subtype, y = expr, fill = subtype)) +
      geom_boxplot() +
      scale_fill_manual(values = lund_colors[subtype_order]) +
      labs(
        title = gene,
        subtitle = paste0("FDR p = ", padj)
      ) +
      theme_minimal() +
      theme(
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
      )
    plots[[i]] <- p
  }
  
  n_col <- ceiling(sqrt(length(plots)))
  n_row <- ceiling(length(plots) / n_col)

  plots[[length(plots) - n_col + 1]] <- plots[[length(plots) - n_col + 1]] +
    labs(x = "Subtype")
  
  plots[[1]] <- plots[[1]] + labs(y = "Expression")
  
  facet_plot = wrap_plots(plots, ncol = n_col) &
    plot_annotation(
      theme = theme(
        plot.margin = margin(5, 5, 5, 5)
      )
    )
  
  #produce volcano plot
  # Prepare data
  df <- anova_df
  
  # Add hub_score and is_signature to df
  metrics <- node_metrics
  df$hub_score <- metrics$hub_score[match(df$gene, metrics$name)]
  df$is_signature <- metrics$is_signature[match(df$gene, metrics$name)]
  
  # Identify top 10 hubs
  top_hub_cutoff <- sort(df$hub_score, decreasing = TRUE)[10]
  df$top_hub <- df$hub_score >= top_hub_cutoff
  
  # Set color: red for top 10 hubs, black otherwise
  df$color <- ifelse(df$top_hub, "red", "black")
  # Set shape: triangle for signature genes, circle otherwise
  df$shape <- ifelse(df$is_signature, 17, 16)
  
  volcano_plot = ggplot(df, aes(x = max_diff, y = -log10(padj), label = gene)) +
    geom_point(aes(color = color, shape = factor(shape)), size = 3) +
    scale_color_identity() +
    scale_shape_manual(values = c("16" = 16, "17" = 17),
                       labels = c("16" = "Other", "17" = "Signature")) +
    geom_text(aes(label = ifelse(padj < 0.05, as.character(gene), "")),
              vjust = 1.5, size = 3, show.legend = FALSE) +
    labs(
      x = "Max difference between subtypes",
      y = "-log10(adjusted p-value)",
      title = "Volcano plot of ANOVA results",
      color = "Top 10 Hub",
      shape = "Signature Gene"
    ) +
    theme_minimal()
  
  #blacklsit genes that are non-significant
  blacklist <- setdiff(extended_genes, results$sig_genes$gene)
  
  results = list(sig_genes = sig_genes,
                 blacklist_genes = blacklist,
                 anova_results = anova_df, 
                 posthoc_results = posthoc_results, 
                 facet_plot = facet_plot,
                 volcano_plot = volcano_plot)
  
  print(results$plot)
  return(results)
}