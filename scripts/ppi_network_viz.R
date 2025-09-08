#load data
load(file = "data/ppi_network/ppi_networks.Rdata")
load(file = "data/hubs/subtype_hubs.Rdata")

#source funciton
source('functions/plot_ppi_network.R')

#generate figures
plot_ppi_network("uro", top_n = 20, top_hubs_highlight = 5)
plot_ppi_network("gu", top_n = 20, top_hubs_highlight = 5)
plot_ppi_network("basq", top_n = 20, top_hubs_highlight = 5)
plot_ppi_network("mes", top_n = 20, top_hubs_highlight = 5)
plot_ppi_network("scne", top_n = 20, top_hubs_highlight = 5)

# Create combined figure with all subtypes
library(gridExtra)
library(grid)

# Generate plots for all subtypes and store them
subtype_plots <- list()

for(subtype in names(ppi_networks)) {
  # Create plot without showing it immediately
  p <- plot_ppi_network(subtype, 
                        top_n = 20, 
                        top_hubs_highlight = 5,
                        show_labels = FALSE)  # Turn off labels for cleaner combined view
  
  # Store the plot
  subtype_plots[[subtype]] <- p
}

# Arrange all plots in a grid
combined_plot <- grid.arrange(
  subtype_plots[["uro"]], 
  subtype_plots[["gu"]], 
  subtype_plots[["basq"]], 
  subtype_plots[["mes"]], 
  subtype_plots[["scne"]], 
  ncol = 3, 
  nrow = 2,
  top = textGrob("PPI Networks by Bladder Cancer Subtype", 
                 gp = gpar(fontsize = 16, fontface = "bold"))
)

# Print the combined plot
print(combined_plot)

# Alternative: Using patchwork for better control
library(patchwork)

# Create combined plot with patchwork - fix for ggraph compatibility
combined_patchwork <- wrap_plots(
  subtype_plots[["uro"]], 
  subtype_plots[["gu"]], 
  subtype_plots[["basq"]], 
  subtype_plots[["mes"]], 
  subtype_plots[["scne"]], 
  ncol = 3, nrow = 2
) + 
plot_annotation(
  title = "PPI Networks by Bladder Cancer Subtype",
  theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5, family = "Arial"))
)

# Print the patchwork version
print(combined_patchwork)

# Save the combined figure as PDF with explicit dimensions
pdf("combined_ppi_networks.pdf", width = 18, height = 12)
print(combined_patchwork)
dev.off()

# Alternative: Save each plot separately then combine
# This is more reliable for ggraph plots
cairo_pdf("combined_ppi_networks_alt.pdf", width = 18, height = 12, family = "Arial")
print(combined_patchwork)
dev.off()

# Backup approach: Manual grid layout if patchwork still fails
if (FALSE) {  # Set to TRUE if needed
  # Create a more robust combined plot using grid
  library(grid)
  library(gridExtra)
  
  pdf("combined_ppi_networks_grid.pdf", width = 18, height = 12)
  
  # Create layout
  grid.newpage()
  
  # Add title
  grid.text("PPI Networks by Bladder Cancer Subtype", 
            x = 0.5, y = 0.95, 
            gp = gpar(fontsize = 16, fontface = "bold"))
  
  # Define plotting areas
  pushViewport(viewport(x = 0.17, y = 0.7, width = 0.3, height = 0.4))
  print(subtype_plots[["uro"]], newpage = FALSE)
  popViewport()
  
  pushViewport(viewport(x = 0.5, y = 0.7, width = 0.3, height = 0.4))
  print(subtype_plots[["gu"]], newpage = FALSE)
  popViewport()
  
  pushViewport(viewport(x = 0.83, y = 0.7, width = 0.3, height = 0.4))
  print(subtype_plots[["basq"]], newpage = FALSE)
  popViewport()
  
  pushViewport(viewport(x = 0.17, y = 0.3, width = 0.3, height = 0.4))
  print(subtype_plots[["mes"]], newpage = FALSE)
  popViewport()
  
  pushViewport(viewport(x = 0.5, y = 0.3, width = 0.3, height = 0.4))
  print(subtype_plots[["scne"]], newpage = FALSE)
  popViewport()
  
  dev.off()
}
