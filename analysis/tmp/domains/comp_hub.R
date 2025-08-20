# Function to plot all subtypes using the same PPI data
plot_all_subtypes <- function(ppi_data, n_hubs = 10, n_genes = 80, save_plots = TRUE) {
  
  plots <- list()
  
  for (subtype in names(expr_subtypes)) {
    cat("\n", "="*50, "\n")
    plots[[subtype]] <- plot_cell_cycle_ppi(
      ppi_data = ppi_data,
      subtype = subtype,
      n_hubs = n_hubs,
      n_genes = n_genes,
      show_labels = FALSE,
      save_plot = save_plots
    )
    print(plots[[subtype]])
  }
  
  return(plots)
}

# Function to quickly compare different hub numbers
compare_hub_numbers <- function(ppi_data, subtype = "uro", hub_counts = c(5, 10, 15)) {
  
  plots <- list()
  
  for (n_hubs in hub_counts) {
    cat("\nPlotting with", n_hubs, "hubs...\n")
    plots[[paste0("hubs_", n_hubs)]] <- plot_cell_cycle_ppi(
      ppi_data = ppi_data,
      subtype = subtype,
      n_hubs = n_hubs,
      n_genes = 80,
      show_labels = FALSE
    )
    print(plots[[paste0("hubs_", n_hubs)]])
  }
  
  return(plots)
}
