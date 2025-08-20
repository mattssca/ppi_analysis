#load packages
library(STRINGdb)
library(igraph)

#source scripts
source("functions/get_ppi.R")

#get data
load("data/bio_markers/subtype_biomarkers.Rdata")
load("data/expr/expr_subtypes.Rdata")

#create subtype-specific networks using the biomarker genes
ppi_networks <- list()

for(subtype_name in names(subtype_biomarkers)) {
  cat("Creating network for", subtype_name, "with", length(subtype_biomarkers[[subtype_name]]), "genes\n")
  
  #get the gene list for this subtype
  genes_for_subtype <- subtype_biomarkers[[subtype_name]]
  
  #subset expression data to only include these genes
  #make sure genes exist in the expression data
  available_genes <- intersect(genes_for_subtype, rownames(expr_subtypes[[subtype_name]]))
  
  if(length(available_genes) > 5) {  # Ensure minimum network size
    #create expression matrix with only the biomarker genes
    expr_subset <- expr_subtypes[[subtype_name]][available_genes, ]
    
    #create PPI network using the subsetted expression data
    ppi_networks[[subtype_name]] <- run_stringdb(expr_data = expr_subset)
    
  } else {
    cat("Warning: Only", length(available_genes), "genes available for", subtype_name, "- skipping\n")
  }
}

#check what networks were created
sapply(ppi_networks, function(x) if(is.null(x)) 0 else vcount(x))

#save results
save(ppi_networks, file = "data/ppi_network/ppi_networks.Rdata")
