#load data
load("data/hubs/hub_matrix.Rdata")
load("data/hubs/subtype_hubs.Rdata")
load("data/hubs/top_hubs_comparison.Rdata")

#load pacakges
library(pheatmap)
library(RColorBrewer)

################HEATMAP HUB Matrix#####################
#visualize hub_matrix as a heatmap
pheatmap(hub_matrix, 
         color = c("white", "red"),
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         main = "Hub Gene Specificity Across Subtypes",
         legend_breaks = c(0, 1),
         legend_labels = c("Not Hub", "Hub"))


################Scatterplot top 10####################
# Scatter plot: expression vs hub importance
# Combine all hub data
all_hubs <- bind_rows(subtype_hubs) %>%
  filter(hub_score <= 10)  # Top 10 for clarity

#plot
ggplot(all_hubs, aes(x = mean_expr, y = hub_score, color = subtype)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = gene), size = 3, max.overlaps = 10) +
  scale_color_manual(values = lund_colors) +
  labs(title = "Hub Centrality vs Expression Level",
       x = "Mean Expression", 
       y = "Hub Score") +
  theme_minimal()

##############Facetted Scatter Plot Subtype###########
# Combine top 30 from each subtype
all_top_20 <- data.frame()

for(subtype_name in names(subtype_hubs)) {
  top_20 <- subtype_hubs[[subtype_name]] %>%
    slice_head(n = 30) %>%
    mutate(subtype = subtype_name)
  
  all_top_20 <- rbind(all_top_20, top_20)
}

# Create faceted plot with Lund colors
ggplot(all_top_20, aes(x = mean_expr, y = hub_score, color = subtype)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_text_repel(aes(label = gene), 
                  size = 2.5, 
                  max.overlaps = 20,
                  box.padding = 0.3) +
  facet_wrap(~toupper(subtype), scales = "free", ncol = 3) +
  scale_color_manual(values = lund_colors) +
  labs(title = "Hub Centrality vs Expression Level by Subtype",
       subtitle = "Top 20 hub genes per subtype",
       x = "Mean Expression", 
       y = "Hub Score") +
  theme_minimal() +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(size = 16, face = "bold"),
        strip.background = element_rect(fill = "gray90", color = "white"))
