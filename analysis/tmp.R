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

#alternatively, remove Mes
subtype_vector_class_7_no_mes <- subtype_vector_class_7[subtype_vector_class_7 != "Mes"]
subtype_vector_class_7_stripped <- subtype_vector_class_7[!subtype_vector_class_7 %in% c("ScNE", "Mes")]
subtype_vector_class_7_stripped <- subtype_vector_class_7_stripped[colnames(expr_z_sub)]
subtype_vector_class_7_stripped <- subtype_vector_class_7_stripped[!is.na(subtype_vector_class_7_stripped)]

#select tat1 samples
these_samples = UROSCANSEQMetadata2025_01_16 %>% 
  rownames_to_column("sample_id") %>% 
  filter(Stage %in% c("Ta", "T1")) %>% 
  pull(sample_id)

samples = names(subtype_vector_class_7_stripped)

#z score transformation of expression values
expr_log2 <- log2(uroscanseq_non_batch_corrected + 1)
expr_z <- t(scale(t(expr_log2)))

expr_z_sub = as.data.frame(expr_z) %>% 
  select(any_of(samples))

#run function
ERBB_nodes = expand_and_plot_signature_network(expr_data = expr_z_sub,
                                  lund_genes = lund_genes,
                                  add_lund_genes = TRUE,
                                  node_size = 25,
                                  plot_width = 20,
                                  plot_height = 20, 
                                  color_scale_global = FALSE,
                                  return_data = TRUE,
                                  expr_summary = "mean",
                                  node_color = "lund",
                                  node_degree = FALSE,
                                  theme = "light",
                                  show_labels = TRUE,
                                  max_added_genes = 80, 
                                  layout_method = "fr",
                                  min_degree = 1, 
                                  string_score_threshold = 500,  
                                  subtype_vector = subtype_vector_class_7_stripped, 
                                  signature_list = lundtax_signatures, 
                                  signature_name = "ERBB",
                                  verbose = TRUE, 
                                  out_dir = "viz/tat1/ERBB/")



















####################################################################################################
#wraper for all signatures
for (sig_name in names(lundtax_signatures)) {
  out_dir <- file.path("viz/networks", sig_name)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  expand_and_plot_signature_network(
    expr_data = expr_mat_uroscanseq,
    max_added_genes = 40,
    layout_method = "fr",
    min_degree = 1,
    string_score_threshold = 500,
    subtype_vector = subtype_vector,
    signature_list = lundtax_signatures,
    signature_name = sig_name,
    verbose = FALSE,
    out_dir = out_dir
  )
}