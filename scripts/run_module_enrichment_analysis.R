##### PART 1 - PPI re-creation with module node status
library(igraph)
library(tidygraph)
library(ggraph)

#subset graph object
g <- erbb_big_network$igraph_object

#detect communities (e.g., Louvain)
comm <- cluster_louvain(g)
modules <- membership(comm)

#create tidygraph object with module infomration
tg <- as_tbl_graph(g) %>%
  mutate(module = as.factor(modules))

#set module colurs
my_colors <- c("#59AC77", "#3B0270", "#E45A92", "#EF7722", "#9A3F3F")


#draw ppi network with module information
p = ggraph(tg, layout = erbb_big_network_layout) +
  geom_edge_link(color = "grey", alpha = 0.5, width = 0.1) +
  geom_node_point(aes(color = module), size = 30) +
  scale_color_manual(values = my_colors) +
  geom_label(aes(x = x, y = y, label = name), fill = "white", color = "black", size = 6, fontface = "bold", label.size = 0.2) +
  theme_void() +
  theme_bw(base_size = 15) +
  labs(title = "ERBB Big Network Modules (Louvain Community Detection)") +
  theme(
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(),
    axis.text.x = element_blank(), 
    panel.grid = element_blank(),
    axis.text.y = element_blank(), 
    axis.ticks = element_blank(),
    plot.title.position = "plot",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )

#export ppi netowork
CairoPDF("analysis_sets/ERBB/figures/big_network/test.pdf", width = 30, height = 30)
print(p)
dev.off()



#### PART 2 - Module enrichment analysis
#load packages
library(clusterProfiler)
library(org.Hs.eg.db)

#for each module, run enrichment
enrich_results <- lapply(unique(modules), function(mod) {
  genes_in_mod <- names(modules)[modules == mod]
  enrichGO(
    gene         = genes_in_mod,
    OrgDb        = org.Hs.eg.db,
    keyType      = "SYMBOL",
    ont          = "BP",
    pAdjustMethod= "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable     = TRUE
  )
})

names(enrich_results) <- paste0("Module_", unique(modules))

#view top GO terms for Module 1
head(as.data.frame(enrich_results[["Module_1"]]))

#dotplot for a module
dotplot(enrich_results[["Module_1"]], showCategory = 10)
dotplot(enrich_results[["Module_2"]], showCategory = 10)
dotplot(enrich_results[["Module_3"]], showCategory = 10)
dotplot(enrich_results[["Module_4"]], showCategory = 10)

