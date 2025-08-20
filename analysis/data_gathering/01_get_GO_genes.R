#load packages
library(msigdbr)
library(dplyr)

#retrieve all human GO Biological Process gene sets
msig_go <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")


#CELL CYCLE
#filter for your two GO terms by their gene set names
cell_cycle_genes <- msig_go %>% 
  filter(gs_name == "GOBP_CELL_CYCLE")

cell_cycle_process_genes <- msig_go %>% 
  filter(gs_name == "GOBP_CELL_CYCLE_PROCESS")

#combine and deduplicate gene symbols
cell_cycle_go <- unique(c(cell_cycle_genes$gene_symbol, cell_cycle_process_genes$gene_symbol))

#save object
save(cell_cycle_go,file = "data/GO/cell_cycle_go.Rdata")

#LUMINAL DIFFERENTIATION
#filter for relevant GO terms (adjust names as needed)
luminal_diff_genes <- msig_go %>% 
  filter(gs_name %in% c("GOBP_EPITHELIAL_CELL_DIFFERENTIATION", 
                        "GOBP_EPITHELIAL_CELL_FATE_COMMITMENT",
                        "GOBP_EPITHELIAL_CELL_DEVELOPMENT"))

#deduplicate gene symbols
luminal_diff_go <- unique(luminal_diff_genes$gene_symbol)

#save object
save(luminal_diff_go, file = "data/GO/luminal_diff_go.Rdata")
