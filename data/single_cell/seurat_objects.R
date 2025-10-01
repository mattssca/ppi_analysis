library(GEOquery)
library(Seurat)

#normal bladder
# Get supplementary files (the actual single-cell data)
getGEOSuppFiles("GSE129845", makeDirectory = TRUE)

# List downloaded files
list.files("GSE129845/", full.names = TRUE)

# Get the actual file paths
all_files <- list.files("GSE129845/", full.names = TRUE)

# Find unique sample identifiers
sample_ids <- unique(gsub(".*_(Sample[0-9]+)_.*", "\\1", basename(all_files)))
sample_ids <- sample_ids[grepl("Sample", sample_ids)]

print("Found samples:")
print(sample_ids)

# Load each sample
seurat_list <- list()

for(sample_id in sample_ids) {
  
  # Find files for this sample
  sample_files <- all_files[grepl(sample_id, all_files)]
  
  matrix_file <- sample_files[grepl("matrix.mtx", sample_files)]
  barcode_file <- sample_files[grepl("barcodes.tsv", sample_files)]
  genes_file <- sample_files[grepl("genes.tsv", sample_files)]
  
  cat("Loading", sample_id, ":\n")
  cat("Matrix:", matrix_file, "\n")
  cat("Barcodes:", barcode_file, "\n")
  cat("Genes:", genes_file, "\n\n")
  
  # Read 10X data
  data <- ReadMtx(
    mtx = matrix_file,
    cells = barcode_file,
    features = genes_file,
    feature.column = 2
  )
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(
    counts = data,
    project = sample_id,
    min.cells = 3,
    min.features = 200
  )
  
  seurat_list[[sample_id]] <- seurat_obj
}

# Merge all samples
bladder_combined <- merge(seurat_list[[1]], 
                          y = seurat_list[2:length(seurat_list)], 
                          add.cell.ids = names(seurat_list))

print(bladder_combined)

bladder_normal = bladder_combined

#bladder cancer
# Download the raw single-cell data files
getGEOSuppFiles("GSE130001", makeDirectory = TRUE)

# Check what files were downloaded
supp_files <- list.files("GSE130001/", full.names = TRUE)
print(supp_files)

# Extract any tar or compressed files
if(any(grepl("\\.tar$", supp_files))) {
  tar_file <- supp_files[grepl("\\.tar$", supp_files)]
  untar(tar_file, exdir = "GSE130001/")
}

# Look for single-cell data files
all_files <- list.files("GSE130001/", full.names = TRUE, recursive = TRUE)
print(all_files)


# First, decompress the .gz files
gz_files <- all_files[grepl("\\.gz$", all_files)]
for(file in gz_files) {
  R.utils::gunzip(file, overwrite = TRUE)
}

# List decompressed files
decompressed_files <- list.files("GSE130001/", full.names = TRUE)
print(decompressed_files)

# Load each sample
sample_names <- c("sample1", "sample2")
seurat_list <- list()

for(i in 1:2) {
  sample_name <- sample_names[i]
  sample_id <- paste0("GSM372917", 7+i)
  
  # File paths
  matrix_file <- paste0("GSE130001/", sample_id, "_", sample_name, "_matrix.mtx")
  barcode_file <- paste0("GSE130001/", sample_id, "_", sample_name, "_barcodes.tsv")
  genes_file <- paste0("GSE130001/", sample_id, "_", sample_name, "_genes.tsv")
  
  cat("Loading", sample_name, ":\n")
  cat("Files:", matrix_file, barcode_file, genes_file, "\n\n")
  
  # Read 10X data
  data <- ReadMtx(
    mtx = matrix_file,
    cells = barcode_file,
    features = genes_file,
    feature.column = 2
  )
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(
    counts = data,
    project = paste("BlCa", sample_name),  # Bladder Cancer
    min.cells = 3,
    min.features = 200
  )
  
  seurat_list[[sample_name]] <- seurat_obj
}

# Merge samples
bladder_cancer_seurat <- merge(seurat_list[["sample1"]], 
                               y = seurat_list[["sample2"]], 
                               add.cell.ids = c("sample1", "sample2"))

print(bladder_cancer_seurat)
cat("Bladder cancer cells:", ncol(bladder_cancer_seurat), "\n")
cat("Genes:", nrow(bladder_cancer_seurat), "\n")

bladder_cancer = bladder_cancer_seurat

save(bladder_cancer, file = "data/single_cell/bladder_cancer_seurat.Rdata")
save(bladder_normal, file = "data/single_cell/bladder_normal_seurat.Rdata")

#standard single-cell processing
bladder_cancer <- NormalizeData(bladder_cancer)
bladder_cancer <- FindVariableFeatures(bladder_cancer)
bladder_cancer <- ScaleData(bladder_cancer)
bladder_cancer <- RunPCA(bladder_cancer)

#number of PCs directly
n_pcs <- ncol(Embeddings(bladder_cancer, "pca"))

#Use the available number of PCs
max_dims <- min(15, n_pcs)

#resume standard single-cell processing
bladder_cancer <- FindNeighbors(bladder_cancer, dims = 1:max_dims)
bladder_cancer <- FindClusters(bladder_cancer, resolution = 0.5)
bladder_cancer <- RunUMAP(bladder_cancer, dims = 1:max_dims)

bladder_cancer_processed = bladder_cancer

#repeat
#standard single-cell processing
bladder_normal <- NormalizeData(bladder_normal)
bladder_normal <- FindVariableFeatures(bladder_normal)
bladder_normal <- ScaleData(bladder_normal)
bladder_normal <- RunPCA(bladder_normal)

#number of PCs directly
n_pcs <- ncol(Embeddings(bladder_normal, "pca"))

#Use the available number of PCs
max_dims <- min(15, n_pcs)

#resume standard single-cell processing
bladder_normal <- FindNeighbors(bladder_normal, dims = 1:max_dims)
bladder_normal <- FindClusters(bladder_normal, resolution = 0.5)
bladder_normal <- RunUMAP(bladder_normal, dims = 1:max_dims)

bladder_normal_processed = bladder_normal

save(bladder_cancer_processed, file = "data/single_cell/bladder_cancer_seurat_processed.Rdata")
save(bladder_normal_processed, file = "data/single_cell/bladder_normal_seurat_processed.Rdata")
