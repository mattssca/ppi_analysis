#load data
load("data/lundtax_signatures.Rdata")
load("data/expr/expr_mat_uroscanseq.Rdata")
load("data/predicted/subtype_vector.Rdata")
load("C:/Users/matts/Desktop/ppi_analysis/data/predicted/predicted_full.Rdata")
load("C:/Users/matts/Desktop/uroscanseq_analysis/data/raw/uroscanseq_non_batch_corrected_754.Rdata")
load("C:/Users/matts/Desktop/uroscanseq_analysis/data/metadata/UROSCANSEQMetadata2025_01_16.Rdata")

#source script
source("analysis/expand_and_plot_signature_network.R")

#create new signature as a combination of late and early genes
combined_genes <- unique(c(lundtax_signatures[["Early_CC"]], lundtax_signatures[["Late_CC"]]))
lundtax_signatures[["Combined_CC"]] <- combined_genes

#select subtype class
subtype_vector_class_7 <- predicted$predictions_7classes

#select tat1 samples
these_samples = UROSCANSEQMetadata2025_01_16 %>% 
  rownames_to_column("sample_id") %>% 
  filter(Stage %in% c("Ta", "T1")) %>% 
  pull(sample_id)

#alternatively, remove Mes
subtype_vector_class_7_stripped <- subtype_vector_class_7[!subtype_vector_class_7 %in% c("ScNE", "Mes")]
subtype_vector_class_7_stripped <- subtype_vector_class_7_stripped[these_samples]
subtype_vector_class_7_stripped <- subtype_vector_class_7_stripped[!is.na(subtype_vector_class_7_stripped)]

#z score transformation of expression values
expr_log2 <- log2(uroscanseq_non_batch_corrected + 1)
expr_z <- t(scale(t(expr_log2)))

expr_z_sub = as.data.frame(expr_z) %>% 
  dplyr::select(any_of(names(subtype_vector_class_7_stripped)))

#run function
erbb_nodes_small = expand_and_plot_signature_network(expr_data = expr_z_sub,
                                                     node_size = 25,
                                                     plot_width = 15,
                                                     plot_height = 15, 
                                                     color_scale_global = FALSE,
                                                     return_data = TRUE,
                                                     expr_summary = "mean",
                                                     node_color = "lund",
                                                     node_degree = FALSE,
                                                     theme = "light",
                                                     show_labels = TRUE,
                                                     max_added_genes = 20, 
                                                     layout_method = "kk",
                                                     min_degree = 1, 
                                                     string_score_threshold = 500,  
                                                     subtype_vector = subtype_vector_class_7_stripped, 
                                                     signature_list = lundtax_signatures, 
                                                     signature_name = "ERBB",
                                                     verbose = TRUE, 
                                                     out_dir = "viz/tat1/ERBB/complete_analysis/")

erbb_nodes_large = expand_and_plot_signature_network(expr_data = expr_z_sub,
                                                     node_size = 25,
                                                     plot_width = 15,
                                                     plot_height = 15, 
                                                     color_scale_global = FALSE,
                                                     return_data = TRUE,
                                                     expr_summary = "mean",
                                                     node_color = "lund",
                                                     node_degree = FALSE,
                                                     theme = "light",
                                                     show_labels = TRUE,
                                                     max_added_genes = 200, 
                                                     layout_method = "kk",
                                                     min_degree = 1, 
                                                     string_score_threshold = 500,  
                                                     subtype_vector = subtype_vector_class_7_stripped, 
                                                     signature_list = lundtax_signatures, 
                                                     signature_name = "ERBB",
                                                     verbose = TRUE, 
                                                     out_dir = "viz/tat1/ERBB/complete_analysis/")
