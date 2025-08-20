# Updated function with correct attributes
get_protein_domains <- function(gene_list, dataset = "hsapiens_gene_ensembl") {
  
  # Connect to Ensembl
  mart <- useMart("ensembl", dataset = dataset)
  
  # Get domain information with correct attribute names
  domains <- getBM(
    attributes = c(
      "hgnc_symbol",
      "interpro_short_description", 
      "interpro_description",
      "interpro",
      "pfam",
      "pfam_start",
      "pfam_end"
    ),
    filters = "hgnc_symbol",
    values = gene_list,
    mart = mart
  )
  
  return(domains)
}

# Alternative approach if Pfam descriptions aren't available
get_protein_domains_simple <- function(gene_list, dataset = "hsapiens_gene_ensembl") {
  
  # Connect to Ensembl
  mart <- useMart("ensembl", dataset = dataset)
  
  # Get domain information - simplified version
  domains <- getBM(
    attributes = c(
      "hgnc_symbol",
      "interpro_short_description", 
      "interpro_description",
      "interpro"
    ),
    filters = "hgnc_symbol",
    values = gene_list,
    mart = mart
  )
  
  return(domains)
}

# Try the updated function
subtype_domains <- list()

for(subtype in names(subtype_hubs)) {
  # Get top 50 hub genes for domain analysis
  hub_genes <- subtype_hubs[[subtype]]$gene[1:50]
  
  cat("Getting domains for", subtype, "subtype...\n")
  
  tryCatch({
    subtype_domains[[subtype]] <- get_protein_domains(hub_genes) %>%
      filter(!is.na(interpro_short_description)) %>%
      distinct()
  }, error = function(e) {
    cat("Error with full function, trying simplified version...\n")
    subtype_domains[[subtype]] <- get_protein_domains_simple(hub_genes) %>%
      filter(!is.na(interpro_short_description)) %>%
      distinct()
  })
}

# Save results
save(subtype_domains, file = "data/domains/subtype_domains.Rdata")
