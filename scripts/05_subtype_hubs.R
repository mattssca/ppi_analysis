#load packages
library(igraph)
library(dplyr)

#source scripts
source("functions/get_subtype_hubs.R")

#load data
load("data/ppi_network/ppi_networks.Rdata")

#run hub analysis on your subtype-specific networks
subtype_hubs <- analyze_subtype_hubs(ppi_networks, top_n = 20)

#save results
save(subtype_hubs, file = "data/hubs/subtype_hubs.Rdata")

#compare top hubs across subtypes
cat("\nTop 5 hubs per subtype:\n")
top_hubs_comparison <- lapply(subtype_hubs, function(x) {
  data.frame(
    gene = x$gene[1:5],
    hub_score = round(x$hub_score[1:5], 2),
    mean_expr = round(x$mean_expr[1:5], 2)
  )
})

print(top_hubs_comparison)

#save top 5 hubs for each subtype
save(top_hubs_comparison, file = "data/hubs/top_hubs_comparison.Rdata")

#check for unique vs shared hubs across subtypes
all_top_hubs <- unique(unlist(lapply(subtype_hubs, function(x) x$gene[1:10])))
hub_matrix <- matrix(0, nrow = length(all_top_hubs), ncol = length(subtype_hubs))
rownames(hub_matrix) <- all_top_hubs
colnames(hub_matrix) <- names(subtype_hubs)

for(subtype in names(subtype_hubs)) {
  top_genes <- subtype_hubs[[subtype]]$gene[1:10]
  hub_matrix[top_genes, subtype] <- 1
}

cat("\nHub specificity (1 = hub in subtype, 0 = not a hub):\n")
print(hub_matrix[rowSums(hub_matrix) <= 2, ])  #show genes that are hubs in ≤2 subtypes

#save results
save(hub_matrix, file = "data/hubs/hub_matrix.Rdata")
