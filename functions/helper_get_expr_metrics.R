# Function to get expression summary statistics
get_expr_summary <- function(expr_data, subtype_name) {
  # Get all expression values
  all_values <- as.numeric(as.matrix(expr_data))
  all_values <- all_values[!is.na(all_values)]
  
  # Calculate statistics
  summary_stats <- data.frame(
    subtype = subtype_name,
    samples = ncol(expr_data),
    genes = nrow(expr_data),
    min = min(all_values),
    q25 = quantile(all_values, 0.25),
    median = median(all_values),
    mean = mean(all_values),
    q75 = quantile(all_values, 0.75),
    max = max(all_values),
    sd = sd(all_values)
  )
  
  return(summary_stats)
}

# Get summaries for each subtype
expr_summaries <- bind_rows(
  get_expr_summary(expr_subtypes$uro, "uro"),
  get_expr_summary(expr_subtypes$gu, "gu"),
  get_expr_summary(expr_subtypes$basq, "basq"),
  get_expr_summary(expr_subtypes$mes, "mes"),
  get_expr_summary(expr_subtypes$scne, "scne")
)

print(expr_summaries)

# Look at gene-level mean expression distributions
gene_means_summary <- data.frame(
  subtype = names(expr_subtypes),
  mean_expr_min = sapply(expr_subtypes, function(x) min(rowMeans(x, na.rm = TRUE))),
  mean_expr_q25 = sapply(expr_subtypes, function(x) quantile(rowMeans(x, na.rm = TRUE), 0.25)),
  mean_expr_median = sapply(expr_subtypes, function(x) median(rowMeans(x, na.rm = TRUE))),
  mean_expr_q75 = sapply(expr_subtypes, function(x) quantile(rowMeans(x, na.rm = TRUE), 0.75)),
  mean_expr_max = sapply(expr_subtypes, function(x) max(rowMeans(x, na.rm = TRUE)))
)

print("Gene-level mean expression ranges:")
print(gene_means_summary)

# Check cell cycle gene expression specifically
cell_cycle_expr_summary <- data.frame(
  subtype = names(expr_subtypes),
  cc_genes_present = sapply(expr_subtypes, function(x) sum(cell_cycle_go %in% rownames(x))),
  cc_mean_expr_median = sapply(expr_subtypes, function(x) {
    cc_genes_in_data <- intersect(cell_cycle_go, rownames(x))
    if(length(cc_genes_in_data) > 0) {
      median(rowMeans(x[cc_genes_in_data, , drop = FALSE], na.rm = TRUE))
    } else NA
  }),
  cc_mean_expr_q75 = sapply(expr_subtypes, function(x) {
    cc_genes_in_data <- intersect(cell_cycle_go, rownames(x))
    if(length(cc_genes_in_data) > 0) {
      quantile(rowMeans(x[cc_genes_in_data, , drop = FALSE], na.rm = TRUE), 0.75)
    } else NA
  })
)

print("Cell cycle gene expression summary:")
print(cell_cycle_expr_summary)

# Plot histograms of mean expression for cell cycle genes
par(mfrow = c(2, 3))
for(subtype in names(expr_subtypes)) {
  cc_genes_in_data <- intersect(cell_cycle_go, rownames(expr_subtypes[[subtype]]))
  if(length(cc_genes_in_data) > 0) {
    mean_expr <- rowMeans(expr_subtypes[[subtype]][cc_genes_in_data, , drop = FALSE], na.rm = TRUE)
    hist(mean_expr, main = paste("Cell Cycle Genes -", subtype), 
         xlab = "Mean Expression", breaks = 20, col = "lightblue")
    abline(v = c(1, 1.5, 2, 3), col = c("red", "orange", "blue", "green"), lty = 2)
  }
}
