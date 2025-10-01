#confidence Level on Confidence Score
#percentile-based approach automatically adapts to each dataset's value ranges...
add_confidence_scoring <- function(classification_df, 
                                   uro_high_detection = 0.3,
                                   uro_high_expression = 0.5,
                                   uro_med_detection = 0.15,
                                   uro_med_expression = 0.2,
                                   nonuro_high_detection = 0.01,
                                   nonuro_high_expression = 0.01,
                                   nonuro_med_detection = 0.05,
                                   nonuro_med_expression = 0.05,
                                   max_cells_threshold = 10000,
                                   max_nonuro_cells = 2000) {
  
  #calculate confidence scores (same as before)
  classification_df$Confidence_Score <- ifelse(
    classification_df$Classification == "Urothelial",
    (classification_df$Detection_Rate * 2) + 
      (pmin(classification_df$Mean_Expression / 2, 1)) + 
      (pmin(classification_df$Cells_Expressing / max_cells_threshold, 1)),
    
    (1 - classification_df$Detection_Rate) * 2 + 
      (1 - pmin(classification_df$Mean_Expression / 0.1, 1)) + 
      (1 - pmin(classification_df$Cells_Expressing / max_nonuro_cells, 1))
  )
  
  #sort by confidence score
  classification_df <- classification_df[order(classification_df$Confidence_Score, decreasing = TRUE), ]
  
  #base confidence level on confidence score percentiles
  #calculate percentiles for each classification type separately
  uro_genes <- classification_df[classification_df$Classification == "Urothelial", ]
  nonuro_genes <- classification_df[classification_df$Classification == "Non-urothelial", ]
  
  #assign confidence levels based on score percentiles
  classification_df$Confidence_Level <- NA
  
  if(nrow(uro_genes) > 0) {
    uro_high_cutoff <- quantile(uro_genes$Confidence_Score, 0.7)  # Top 30%
    uro_med_cutoff <- quantile(uro_genes$Confidence_Score, 0.3)   # Middle 40%
    
    classification_df[classification_df$Classification == "Urothelial", "Confidence_Level"] <- 
      ifelse(uro_genes$Confidence_Score >= uro_high_cutoff, "High",
             ifelse(uro_genes$Confidence_Score >= uro_med_cutoff, "Medium", "Low"))
  }
  
  if(nrow(nonuro_genes) > 0) {
    nonuro_high_cutoff <- quantile(nonuro_genes$Confidence_Score, 0.7)  # Top 30%
    nonuro_med_cutoff <- quantile(nonuro_genes$Confidence_Score, 0.3)   # Middle 40%
    
    classification_df[classification_df$Classification == "Non-urothelial", "Confidence_Level"] <- 
      ifelse(nonuro_genes$Confidence_Score >= nonuro_high_cutoff, "High",
             ifelse(nonuro_genes$Confidence_Score >= nonuro_med_cutoff, "Medium", "Low"))
  }
  
  #print summary
  cat("\n=== CONFIDENCE SUMMARY (Score-Based) ===\n")
  print(table(classification_df$Classification, classification_df$Confidence_Level))
  
  #show correlation
  cat("\nConfidence Score vs Level correlation:\n")
  print(table(classification_df$Confidence_Level, cut(classification_df$Confidence_Score, breaks = 3, labels = c("Low_Score", "Med_Score", "High_Score"))))
  
  return(classification_df)
}

#apply the function
normal_fixed <- add_confidence_scoring(normal_bladder_urothelial, max_cells_threshold = 2000)
cancer_fixed <- add_confidence_scoring(cancer_bladder_urothelial,max_cells_threshold = 1000)

cancer_bladder_urothelial = cancer_fixed
normal_bladder_urothelial = normal_fixed

save(normal_bladder_urothelial, file = "data/single_cell/normal_bladder_urothelial.Rdata")
save(cancer_bladder_urothelial, file = "data/single_cell/cancer_bladder_urothelial.Rdata")
