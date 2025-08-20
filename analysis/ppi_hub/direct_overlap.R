# Get the classifier gene list
classifier_genes <- LundTax2023Classifier::gene_list$hgnc_symbol

# Check overlap for each subtype
hub_classifier_overlap <- list()

for(subtype in names(subtype_hubs)) {
  hub_genes <- subtype_hubs[[subtype]]$gene
  
  # Find overlapping genes
  overlap <- intersect(hub_genes, classifier_genes)
  
  # Calculate enrichment statistics
  hub_classifier_overlap[[subtype]] <- list(
    overlap_genes = overlap,
    n_overlap = length(overlap),
    n_hub_genes = length(hub_genes),
    n_classifier_genes = length(classifier_genes),
    overlap_percent_of_hubs = round(length(overlap) / length(hub_genes) * 100, 1),
    overlap_percent_of_classifier = round(length(overlap) / length(classifier_genes) * 100, 1)
  )
}

# Print summary
for(subtype in names(hub_classifier_overlap)) {
  cat("\n", toupper(subtype), "Subtype:\n")
  cat("Overlapping genes:", paste(hub_classifier_overlap[[subtype]]$overlap_genes, collapse = ", "), "\n")
  cat("Overlap:", hub_classifier_overlap[[subtype]]$n_overlap, "/", hub_classifier_overlap[[subtype]]$n_hub_genes, 
      "hub genes (", hub_classifier_overlap[[subtype]]$overlap_percent_of_hubs, "%)\n")
}

# Fisher's exact test for enrichment
library(broom)

fisher_results <- data.frame()

for(subtype in names(subtype_hubs)) {
  hub_genes <- subtype_hubs[[subtype]]$gene
  overlap <- intersect(hub_genes, classifier_genes)
  
  # Create contingency table
  # Assuming you have a background of all cell cycle genes
  load("data/GO/cell_cycle_go.Rdata")
  
  contingency_table <- matrix(c(
    length(overlap),  # Hub genes that are also classifier genes
    length(hub_genes) - length(overlap),  # Hub genes not in classifier
    length(intersect(classifier_genes, cell_cycle_go)) - length(overlap),  # Classifier genes not hubs
    length(cell_cycle_go) - length(hub_genes) - (length(intersect(classifier_genes, cell_cycle_go)) - length(overlap))  # Neither
  ), nrow = 2)
  
  fisher_test <- fisher.test(contingency_table)
  
  fisher_results <- rbind(fisher_results, data.frame(
    subtype = subtype,
    p_value = fisher_test$p.value,
    odds_ratio = fisher_test$estimate,
    conf_low = fisher_test$conf.int[1],
    conf_high = fisher_test$conf.int[2],
    n_overlap = length(overlap)
  ))
}

print(fisher_results)

library(ggplot2)
library(dplyr)

# Create overlap matrix
overlap_matrix <- matrix(0, nrow = length(classifier_genes), ncol = length(names(subtype_hubs)))
rownames(overlap_matrix) <- classifier_genes
colnames(overlap_matrix) <- names(subtype_hubs)

for(subtype in names(subtype_hubs)) {
  hub_genes <- subtype_hubs[[subtype]]$gene
  overlap_genes <- intersect(hub_genes, classifier_genes)
  overlap_matrix[overlap_genes, subtype] <- 1
}

# Convert to long format for ggplot
overlap_df <- as.data.frame(overlap_matrix) %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "subtype", values_to = "is_hub") %>%
  filter(is_hub == 1)  # Only show overlapping genes

# Plot
ggplot(overlap_df, aes(x = subtype, y = gene, fill = factor(is_hub))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("1" = "red"), name = "Hub Gene") +
  labs(title = "Classifier Genes That Are Also Hub Genes",
       x = "Subtype", y = "Classifier Gene") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1))

library(UpSetR)

# Create list of hub genes per subtype that overlap with classifier
hub_classifier_lists <- list()
for(subtype in names(subtype_hubs)) {
  hub_genes <- subtype_hubs[[subtype]]$gene
  hub_classifier_lists[[subtype]] <- intersect(hub_genes, classifier_genes)
}

# Remove empty lists
hub_classifier_lists <- hub_classifier_lists[sapply(hub_classifier_lists, length) > 0]

if(length(hub_classifier_lists) > 1) {
  upset(fromList(hub_classifier_lists), 
        order.by = "freq",
        main.bar.color = "steelblue",
        sets.bar.color = "darkred",
        text.scale = 1.2,
        mainbar.y.label = "Classifier Gene Intersections",
        sets.x.label = "Classifier Genes per Subtype")
}

# Analyze hub scores of classifier genes vs non-classifier genes
classifier_hub_analysis <- data.frame()

for(subtype in names(subtype_hubs)) {
  subtype_data <- subtype_hubs[[subtype]] %>%
    mutate(
      is_classifier = gene %in% classifier_genes,
      classifier_status = ifelse(is_classifier, "Classifier Gene", "Other Hub Gene")
    )
  
  classifier_hub_analysis <- rbind(classifier_hub_analysis, subtype_data)
}

# Compare hub scores
ggplot(classifier_hub_analysis, aes(x = classifier_status, y = hub_score, fill = classifier_status)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.6) +
  facet_wrap(~subtype) +
  labs(title = "Hub Scores: Classifier vs Non-Classifier Genes",
       y = "Hub Score (lower = more central)",
       x = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Statistical test
library(broom)
hub_score_tests <- classifier_hub_analysis %>%
  group_by(subtype) %>%
  do(test = t.test(hub_score ~ is_classifier, data = .)) %>%
  mutate(p_value = map_dbl(test, ~ .x$p.value),
         mean_diff = map_dbl(test, ~ diff(.x$estimate)))

print(hub_score_tests)

# Compare expression levels of classifier vs non-classifier hub genes
ggplot(classifier_hub_analysis, aes(x = classifier_status, y = mean_expr, fill = classifier_status)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.6) +
  facet_wrap(~subtype) +
  labs(title = "Expression Levels: Classifier vs Non-Classifier Hub Genes",
       y = "Mean Expression",
       x = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Compare different centrality measures
centrality_comparison <- classifier_hub_analysis %>%
  dplyr::select(gene, subtype, is_classifier, degree, betweenness, eigenvector) %>%
  pivot_longer(cols = c(degree, betweenness, eigenvector), 
               names_to = "centrality_measure", 
               values_to = "value")

ggplot(centrality_comparison, aes(x = is_classifier, y = value, fill = is_classifier)) +
  geom_boxplot() +
  facet_grid(centrality_measure ~ subtype, scales = "free_y") +
  labs(title = "Network Centrality: Classifier vs Non-Classifier Genes",
       x = "Is Classifier Gene",
       y = "Centrality Value") +
  theme_minimal()

# Create comprehensive summary
overlap_summary <- data.frame(
  subtype = names(hub_classifier_overlap),
  n_overlap = sapply(hub_classifier_overlap, function(x) x$n_overlap),
  percent_hubs_are_classifier = sapply(hub_classifier_overlap, function(x) x$overlap_percent_of_hubs),
  overlapping_genes = sapply(hub_classifier_overlap, function(x) paste(x$overlap_genes, collapse = ", "))
)

print("Hub-Classifier Overlap Summary:")
print(overlap_summary)

# Overall statistics
total_unique_overlaps <- unique(unlist(sapply(hub_classifier_overlap, function(x) x$overlap_genes)))
cat("\nTotal unique classifier genes found as hubs:", length(total_unique_overlaps))
cat("\nOut of", length(classifier_genes), "total classifier genes")
cat("\nOverall overlap percentage:", round(length(total_unique_overlaps)/length(classifier_genes)*100, 1), "%")
