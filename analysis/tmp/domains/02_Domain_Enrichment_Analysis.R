# Fixed function without Pfam references
analyze_domain_enrichment <- function(subtype_domains) {
  
  domain_analysis <- list()
  
  for(subtype in names(subtype_domains)) {
    
    # Count InterPro domains
    domain_counts <- subtype_domains[[subtype]] %>%
      count(interpro_short_description, sort = TRUE) %>%
      filter(n >= 2)  # Only domains in 2+ proteins
    
    # Count InterPro IDs (more specific)
    interpro_id_counts <- subtype_domains[[subtype]] %>%
      filter(!is.na(interpro)) %>%
      count(interpro, sort = TRUE) %>%
      filter(n >= 2)
    
    domain_analysis[[subtype]] <- list(
      interpro_domains = domain_counts,
      interpro_ids = interpro_id_counts,
      total_proteins = length(unique(subtype_domains[[subtype]]$hgnc_symbol)),
      total_domains = nrow(domain_counts)
    )
  }
  
  return(domain_analysis)
}

# Run the fixed function
domain_enrichment <- analyze_domain_enrichment(subtype_domains)

# Print summary
for(subtype in names(domain_enrichment)) {
  cat("\n", toupper(subtype), "Subtype Domain Summary:\n")
  cat("Proteins analyzed:", domain_enrichment[[subtype]]$total_proteins, "\n")
  cat("Unique domains found:", domain_enrichment[[subtype]]$total_domains, "\n")
  cat("Top domains:\n")
  print(head(domain_enrichment[[subtype]]$interpro_domains, 5))
}
