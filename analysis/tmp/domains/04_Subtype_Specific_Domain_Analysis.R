# Find domains unique to each subtype
find_subtype_specific_domains <- function(domain_enrichment) {
  
  # Get all domains per subtype
  subtype_domain_lists <- list()
  
  for(subtype in names(domain_enrichment)) {
    subtype_domain_lists[[subtype]] <- domain_enrichment[[subtype]]$interpro_domains$interpro_short_description
  }
  
  # Find unique domains
  unique_domains <- list()
  
  for(subtype in names(subtype_domain_lists)) {
    other_domains <- unlist(subtype_domain_lists[names(subtype_domain_lists) != subtype])
    
    unique_domains[[subtype]] <- setdiff(subtype_domain_lists[[subtype]], other_domains)
  }
  
  return(unique_domains)
}

subtype_specific_domains <- find_subtype_specific_domains(domain_enrichment)

# Print unique domains
for(subtype in names(subtype_specific_domains)) {
  cat("\n", toupper(subtype), "- Unique Domains:\n")
  if(length(subtype_specific_domains[[subtype]]) > 0) {
    cat(paste(subtype_specific_domains[[subtype]], collapse = "\n"), "\n")
  } else {
    cat("No unique domains found\n")
  }
}
