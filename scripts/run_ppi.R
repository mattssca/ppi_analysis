run_ppi_analysis = function(my_gene_set = NULL,
                            my_expressions = NULL,
                            my_predictions = NULL,
                            bm_min_expr = 1.5,
                            bm_fold_change = 1,
                            bm_top_n = 200,
                            hub_top = 20,
                            gene_set_name = "My Geneset",
                            plot_ppi_top = 20,
                            plot_ppi_top_hubs_highlight = 5,
                            plot_ppi_labels = TRUE,
                            plot_full_ppi_top = 50,
                            plot_full_ppi_hubs = 10,
                            plot_full_ppi_labels = TRUE,
                            plot_layout = "fr",
                            export_data = FALSE,
                            out_path = NULL){
  
  #load packages
  library(dplyr)
  library(tibble)
  library(pheatmap)
  library(RColorBrewer)
  library(gridExtra)
  library(grid)
  library(patchwork)
  library(igraph)
  library(ggraph)
  library(tidygraph)
  library(ggplot2)
  library(tibble)
  
  
  # Ensure gene_set_name is a single string
  gene_set_name <- as.character(gene_set_name)
  if (length(gene_set_name) > 1) gene_set_name <- paste(gene_set_name, collapse = ", ")

  #source scripts
  source("functions/get_subtype_biomarkers.R")
  source("functions/get_ppi.R")
  source("functions/get_subtype_hubs.R")
  source('functions/plot_ppi_network.R')
  source("functions/plot_full_ppi_plot.R")
  
  
  #define the Lund colors
  lund_colors <- c(
    "uro" = "#3cb44b",
    "gu" = "#4363d8", 
    "basq" = "#CD2626",
    "mes" = "#f58231",
    "scne" = "#A020F0"
  )
  
  
  ## Part 1 - Data Subset
  #subset gene expression data to genes of interest
  expr = my_expressions %>% 
    rownames_to_column("genes") %>% 
    filter(genes %in% my_gene_set) %>% 
    column_to_rownames("genes")

  #get subtype information
  # Ensure predictions_5classes exists
  if(!"predictions_5classes" %in% names(my_predictions)) {
    stop("my_predictions must contain 'predictions_5classes'!")
  }
  these_subtypes = as.data.frame(my_predictions$predictions_5classes) %>% 
    rownames_to_column("sample_id")
  # Rename the column to 'subtype' (assume only one column after rownames)
  colnames(these_subtypes)[2] <- "subtype"
  
  #get sample sets
  uro_samples = these_subtypes %>% 
    filter(subtype == "Uro") %>% 
    pull(sample_id)
  
  gu_samples = these_subtypes %>% 
    filter(subtype == "GU") %>% 
    pull(sample_id)
  
  basq_samples = these_subtypes %>% 
    filter(subtype == "BaSq") %>% 
    pull(sample_id)
  
  mes_samples = these_subtypes %>% 
    filter(subtype == "Mes") %>% 
    pull(sample_id)
  
  scne_samples = these_subtypes %>% 
    filter(subtype == "ScNE") %>% 
    pull(sample_id)
  
  subtypes_samples = list(
    uro = uro_samples,
    gu = gu_samples,
    basq = basq_samples,
    mes = mes_samples,
    scne = scne_samples
  )
  
  #subset to sample groups
  expr_uro = expr %>% dplyr::select(any_of(uro_samples))
  expr_gu = expr %>% dplyr::select(any_of(gu_samples))
  expr_basq = expr %>% dplyr::select(any_of(basq_samples))
  expr_mes = expr %>% dplyr::select(any_of(mes_samples))
  expr_scne = expr %>% dplyr::select(any_of(scne_samples))

  expr_subtypes = list(
    uro = expr_uro,
    gu = expr_gu,
    basq = expr_basq,
    mes = expr_mes,
    scne = expr_scne
  )
  
  ## Part 2 - Gene filter per subtype
  
  #run get_biomarker_genes
  subtype_biomarkers <- list()
  for(subtype_name in names(expr_subtypes)) {
    other_subtypes <- expr_subtypes[names(expr_subtypes) != subtype_name]
    # Only run if there are samples in the subtype
    if(ncol(expr_subtypes[[subtype_name]]) > 0) {
      subtype_biomarkers[[subtype_name]] <- get_biomarker_genes(
        expr_target = expr_subtypes[[subtype_name]],
        expr_others_list = other_subtypes,
        gene_list = my_gene_set,
        min_expr = bm_min_expr,
        fold_change = bm_fold_change,
        top_n = bm_top_n
      )
    } else {
      subtype_biomarkers[[subtype_name]] <- character(0)
      cat("Warning: No samples for subtype", subtype_name, "\n")
    }
  }

  #check results
  sapply(subtype_biomarkers, length)
  
  
  ## Part 3 - Run PPI network analysis
  #create subtype-specific networks using the biomarker genes
  ppi_networks <- list()
  
  for(subtype_name in names(subtype_biomarkers)) {
    cat("Creating network for", subtype_name, "with", length(subtype_biomarkers[[subtype_name]]), "genes\n")
    #get the gene list for this subtype
    genes_for_subtype <- subtype_biomarkers[[subtype_name]]
    #subset expression data to only include these genes
    available_genes <- intersect(genes_for_subtype, rownames(expr_subtypes[[subtype_name]]))
    if(length(available_genes) > 5) {  # Ensure minimum network size
      expr_subset <- expr_subtypes[[subtype_name]][available_genes, , drop=FALSE]
      #create PPI network using the subsetted expression data
      ppi_networks[[subtype_name]] <- run_stringdb(expr_data = expr_subset)
    } else {
      cat("Warning: Only", length(available_genes), "genes available for", subtype_name, "- skipping\n")
      ppi_networks[[subtype_name]] <- NULL
    }
  }
  
  #check what networks were created
  sapply(ppi_networks, function(x) if(is.null(x)) 0 else vcount(x))
  
  
  ## Part 4 - Subtype Hub Identification
  #run hub analysis on your subtype-specific networks
  subtype_hubs <- analyze_subtype_hubs(ppi_networks, top_n = hub_top)
  
  #compare top hubs across subtypes
  cat("\nTop 5 hubs per subtype:\n")
  top_hubs_comparison <- lapply(subtype_hubs, function(x) {
    data.frame(
      gene = x$gene[1:5],
      hub_score = round(x$hub_score[1:5], 2),
      mean_expr = round(x$mean_expr[1:5], 2)
    )
  })
  
  print(top_hubs_comparison)
  
  #check for unique vs shared hubs across subtypes
  all_top_hubs <- unique(unlist(lapply(subtype_hubs, function(x) x$gene[1:10])))
  hub_matrix <- matrix(0, nrow = length(all_top_hubs), ncol = length(subtype_hubs))
  rownames(hub_matrix) <- all_top_hubs
  colnames(hub_matrix) <- names(subtype_hubs)
  
  for(subtype in names(subtype_hubs)) {
    top_genes <- subtype_hubs[[subtype]]$gene[1:min(10, length(subtype_hubs[[subtype]]$gene))]
    hub_matrix[top_genes, subtype] <- 1
  }
  
  cat("\nHub specificity (1 = hub in subtype, 0 = not a hub):\n")
  print(hub_matrix[rowSums(hub_matrix) <= 2, ])
  
  ## Part 5 - Visualization
  ################HEATMAP HUB Matrix#####################
  #visualize hub_matrix as a heatmap
  p_hm_hub_comp = pheatmap(hub_matrix, 
                           color = c("white", "red"),
                           cluster_rows = TRUE,
                           cluster_cols = FALSE,
                           main = paste0("Hub Gene Specificity Across Subtypes (", gene_set_name, ")"),
                           legend_breaks = c(0, 1),
                           legend_labels = c("Not Hub", "Hub"))
  
  ################Scatterplot top 10####################
  #scatter plot: expression vs hub importance
  #combine all hub data
  all_hubs <- bind_rows(subtype_hubs) %>%
    filter(hub_score <= 10)
  
  #plot
  p_scatter_all_hubs = ggplot(all_hubs, aes(x = mean_expr, y = hub_score, color = subtype)) +
                              geom_point(size = 3) +
                              geom_text_repel(aes(label = gene), size = 3, max.overlaps = 10) +
                              scale_color_manual(values = lund_colors) +
                              labs(title = paste0("Hub Centrality vs Expression Level (",my_gene_set, ")"),
                                   x = "Mean Expression", 
                                   y = "Hub Score") +
                              theme_minimal()
  
  ##############Facetted Scatter Plot Subtype###########
  #combine top 30 from each subtype
  all_top_30 <- data.frame()
  
  for(subtype_name in names(subtype_hubs)) {
    top_30 <- subtype_hubs[[subtype_name]] %>%
      slice_head(n = 30) %>%
      mutate(subtype = subtype_name)
    all_top_30 <- rbind(all_top_30, top_30)
  }
  
  #create faceted plot with Lund colors
  p_scatter_hub_facet = ggplot(all_top_30, aes(x = mean_expr, y = hub_score, color = subtype)) +
    geom_point(size = 2.5, alpha = 0.8) +
    geom_text_repel(aes(label = gene), 
                    size = 2.5, 
                    max.overlaps = 20,
                    box.padding = 0.3) +
    facet_wrap(~toupper(subtype), scales = "free", ncol = 3) +
    scale_color_manual(values = lund_colors) +
    labs(title = paste0("Hub Centrality vs Expression Level (",my_gene_set, ")"),
         subtitle = "Top 30 hub genes per subtype",
         x = "Mean Expression", 
         y = "Hub Score") +
    theme_minimal() +
    theme(legend.position = "none",
          strip.text = element_text(face = "bold", size = 12),
          plot.title = element_text(size = 16, face = "bold"),
          strip.background = element_rect(fill = "gray90", color = "white"))
  
  ##################Plot PPI Networks###################
  #generate plots for all subtypes and store them
  subtype_plots <- list()

  # Make ppi_networks available globally for plot_ppi_network
  assign("ppi_networks", ppi_networks, envir = .GlobalEnv)
  assign("subtype_hubs", subtype_hubs, envir = .GlobalEnv)

  for(subtype in names(ppi_networks)) {
    if(!is.null(ppi_networks[[subtype]])) {
      p <- plot_ppi_network(subtype, 
                            top_n = plot_ppi_top, 
                            top_hubs_highlight = plot_ppi_top_hubs_highlight,
                            show_labels = plot_ppi_labels,
                            layout_type = plot_layout)
      subtype_plots[[subtype]] <- p
    }
  }
  
  #create combined plot with patchwork
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
  
  ##Part 6 - Create PPI non-subtype-dependent
  #calculate mean expression for each gene in each subtype
  gene_subtype_expression <- data.frame()
  
  for(subtype in names(expr_subtypes)) {
    # Calculate mean expression for each gene in this subtype
    mean_expr <- rowMeans(expr_subtypes[[subtype]], na.rm = TRUE)
    
    temp_df <- data.frame(
      gene = names(mean_expr),
      subtype = subtype,
      mean_expression = mean_expr
    )
    
    gene_subtype_expression <- rbind(gene_subtype_expression, temp_df)
  }
  
  #for each gene, find which subtype has the highest expression
  gene_dominant_subtype <- gene_subtype_expression %>%
    group_by(gene) %>%
    slice_max(mean_expression, n = 1) %>%
    ungroup() %>%
    select(gene, dominant_subtype = subtype, max_expression = mean_expression)
  
  #check the results
  cat("Genes per dominant subtype:\n")
  table(gene_dominant_subtype$dominant_subtype)
  
  #create PPI network using the full cell cycle expression data
  cat("\nCreating PPI network from expression data...\n")
  full_ppi_network <- run_stringdb(expr_data = expr)
  
  cat("Network created with", vcount(full_ppi_network), "nodes and", ecount(full_ppi_network), "edges\n")
  
  #add subtype dominance information to the network
  #get node names from the network
  network_genes <- V(full_ppi_network)$name
  
  #match with our dominant subtype data
  gene_colors <- gene_dominant_subtype[match(network_genes, gene_dominant_subtype$gene), ]
  
  #add attributes to the network
  V(full_ppi_network)$dominant_subtype <- gene_colors$dominant_subtype
  V(full_ppi_network)$max_expression <- gene_colors$max_expression
  V(full_ppi_network)$node_color <- lund_colors[gene_colors$dominant_subtype]
  
  #handle any missing values (genes not found in subtype data)
  missing_genes <- is.na(V(full_ppi_network)$dominant_subtype)
  if(any(missing_genes)) {
    V(full_ppi_network)$dominant_subtype[missing_genes] <- "unknown"
    V(full_ppi_network)$node_color[missing_genes] <- "gray"
    cat("Warning:", sum(missing_genes), "genes not found in subtype data, colored gray\n")
  }
  
  #create the plots
  full_plot_clean <- create_full_ppi_plot(top_n = plot_full_ppi_top, 
                                          show_labels = plot_full_ppi_labels, 
                                          highlight_top_hubs = plot_full_ppi_hubs,
                                          layout_type = plot_layout)
  
  
  ##Part 7 - Export Plots and Data
  #figures
  # Ensure viz directory exists
  if(!dir.exists("viz")) dir.create("viz")
  cairo_pdf(file.path("viz", paste0(gene_set_name, "_hub_comparison_heatmap.pdf")), width = 7, height = 10, pointsize = 12)
  print(p_hm_hub_comp)
  dev.off()

  cairo_pdf(file.path("viz", paste0(gene_set_name, "_hub_all_subs_scatter.pdf")), width = 10, height = 10, pointsize = 12)
  print(p_scatter_all_hubs)
  dev.off()

  cairo_pdf(file.path("viz", paste0(gene_set_name, "_hub_facet_subs_scatter.pdf")), width = 10, height = 10, pointsize = 12)
  print(p_scatter_hub_facet)
  dev.off()

  cairo_pdf(file.path("viz", paste0(gene_set_name, "_ppi_network_sub_specific_combined.pdf")), width = 20, height = 10, pointsize = 12)
  print(combined_patchwork)
  dev.off()

  cairo_pdf(file.path("viz", paste0(gene_set_name, "_ppi_network_complete.pdf")), width = 20, height = 10, pointsize = 12)
  print(full_plot_clean)
  dev.off()
  
  #data
  if(export_data){
    my_data = list(gene_set = my_gene_set,
                   expr = expr,
                   subtypes_samples = subtypes_samples,
                   expr_subtypes = expr_subtypes,
                   subtype_biomarkers = subtype_biomarkers,
                   ppi_networks = ppi_networks,
                   subtype_hubs = subtype_hubs,
                   top_hubs_comparison = top_hubs_comparison,
                   hub_matrix = hub_matrix,
                   gene_dominant_subtype = gene_dominant_subtype,
                   full_ppi_network = full_ppi_network)
    
    message("Success!")
    return(my_data)
  }
  
  message("Success!")
  
}
