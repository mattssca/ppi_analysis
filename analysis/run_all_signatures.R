# Script to run expand_and_plot_signature_network for all signatures in lundtax_signatures

# Load required libraries and data
# ...existing code to load expr_mat_uroscanseq, subtype_vector, lundtax_signatures...

for (sig_name in names(lundtax_signatures)) {
  out_dir <- file.path("viz/networks", sig_name)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  expand_and_plot_signature_network(
    expr_data = expr_mat_uroscanseq,
    max_added_genes = 40,
    layout_method = "kk",
    min_degree = 1,
    string_score_threshold = 500,
    subtype_vector = subtype_vector,
    signature_list = lundtax_signatures,
    signature_name = sig_name,
    verbose = FALSE,
    out_dir = out_dir
  )
}
