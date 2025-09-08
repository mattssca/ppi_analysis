#load data
load("../uroscanseq_analysis/data/raw/uroscanseq_batch_corrected_572.Rdata")
load("data/lundtax/lundtax_signatures.Rdata")

#zscore transform
expr_full_batchadjusted_zscores = t(scale(t(uroscanseq_batch_corrected)))

#convert back to data frame
expr_full_batchadjusted_zscores = as.data.frame(expr_full_batchadjusted_zscores)

#run predictor
predicted_batchadjusted = lundtax_predict_sub(this_data = expr_full_batchadjusted_zscores, 
                                              log_transform = FALSE, 
                                              adjust = TRUE, 
                                              impute = TRUE, 
                                              include_data = TRUE)

#create subtype vector, excluding Mes and ScNE
subtype_vector = predicted_batchadjusted$predictions_7classes
subtype_vector <- subtype_vector[!subtype_vector %in% c("ScNE", "Mes")]
subtype_vector <- subtype_vector[!is.na(subtype_vector)]

expr_sub_full_batchadjusted_zscores = as.data.frame(expr_full_batchadjusted_zscores) %>% 
  dplyr::select(any_of(names(subtype_vector)))

#rerun the predictor for sanity purposes (should be no Mes or ScNE)
predicted_sub_batchadjusted = lundtax_predict_sub(this_data = expr_sub_full_batchadjusted_zscores, 
                                                  log_transform = FALSE, 
                                                  adjust = TRUE, 
                                                  impute = TRUE, 
                                                  include_data = TRUE)

#create bundled data object for downstream analysis
uroscanseq_data = list(
  expr_df = predicted_sub_batchadjusted$data,
  subtype_7_vector = predicted_sub_batchadjusted$predictions_7classes,
  subtype_5_vector = predicted_sub_batchadjusted$predictions_5classes,
  signature_genes = lundtax_signatures
)

#save analysis object
save(uroscanseq_data, file = "analysis_sets/uroscanseq_data.Rdata")
