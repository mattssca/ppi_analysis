#source function
source("functions/get_subtype_biomarkers.R")

#load data
load("data/expr/expr_subtypes.Rdata")
load("data/GO/cell_cycle_go.Rdata")

#initate an empty list
subtype_biomarkers <- list()

#parameter motivation:
# - min_expr: 1.5 (below Q25, captures moderately expressed genes)
# - fold_change: 1 
# No expression filtering - includes all genes above min_expr threshold
# Network topology will do the filtering - poorly connected genes won't become hubs
# Hub analysis is the real filter - you're only keeping top 20 anyway
# Maximizes gene pool - gives STRING database more genes to work with for better networks
# Equal treatment - all subtypes get similar numbers of genes to work with
# - top_n: 200 (generous but not excessive)

#run function
for(subtype_name in names(expr_subtypes)) {
  other_subtypes <- expr_subtypes[names(expr_subtypes) != subtype_name]
  
  subtype_biomarkers[[subtype_name]] <- get_biomarker_genes(
    expr_target = expr_subtypes[[subtype_name]],
    expr_others_list = other_subtypes,
    gene_list = cell_cycle_go,
    min_expr = 1.5,
    fold_change = 1,
    top_n = 200
  )
}

#check results
sapply(subtype_biomarkers, length)

#export results
save(subtype_biomarkers, file = "data/bio_markers/subtype_biomarkers.Rdata")
