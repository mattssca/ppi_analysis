# Updated plot function
plot_domain_comparison <- function(domain_enrichment) {
  
  lund_colors = LundTax2023Classifier::lund_colors$lund_colors
  
  # Combine all domains
  all_domains <- data.frame()
  
  for(subtype in names(domain_enrichment)) {
    temp_df <- domain_enrichment[[subtype]]$interpro_domains %>%
      slice_head(n = 10) %>%
      mutate(subtype = subtype)
    
    all_domains <- rbind(all_domains, temp_df)
  }
  
  # Plot
  ggplot(all_domains, aes(x = reorder(interpro_short_description, n), 
                          y = n, fill = subtype)) +
    geom_col() +
    facet_wrap(~subtype, scales = "free") +
    coord_flip() +
    scale_fill_manual(values = lund_colors) +
    labs(title = "Top Protein Domains by Subtype",
         x = "Protein Domain",
         y = "Number of Hub Proteins") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8),
          legend.position = "none")
}

plot_domain_comparison(domain_enrichment)
