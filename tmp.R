#load data
load("~/BIOINFORMATICS/git_repos/ppi_analysis/data/lundtax_signatures.Rdata")
load("~/BIOINFORMATICS/git_repos/ppi_analysis/data/expr/expr_mat_uroscanseq.Rdata")
load("~/BIOINFORMATICS/git_repos/ppi_analysis/data/predicted/subtype_vector.Rdata")

#source script
source("analysis/expand_and_plot_signature_network.R")

class_7 = as.data.frame(predicted$predictions_7classes)
subtype_vector_class_7 <- setNames(predicted$predictions_7classes, rownames(predicted$predictions_7classes))

#run function
expand_and_plot_signature_network(expr_data = expr_mat_uroscanseq,
                                  show_labels = FALSE,
                                  lund_colors = FALSE,
                                  max_added_genes = 100, 
                                  layout_method = "kk",
                                  min_degree = 1, 
                                  string_score_threshold = 500,  
                                  subtype_vector = subtype_vector_class_7, 
                                  signature_list = lundtax_signatures, 
                                  signature_name = "UroDiff",
                                  verbose = TRUE, 
                                  out_dir = "viz/networks/test/")

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

