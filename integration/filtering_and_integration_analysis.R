"""
Author: Saranya Canchi
Date Created: 2024-08
Description:
- Performs cell-level filtering, normalization, and integration of CellBender-corrected snRNA-seq data.
- Generates UMAP plots and clustering results to visualize sample integration and cell population structure.
"""

# load libraries
library(Seurat)
library(SeuratWrappers)
library(cluster)
library(ggplot2)
library(future)
options(future.globals.maxSize = 50 * 1024^3)

# set num_cores
num_cores <- future::availableCores()
plan("multisession", workers = num_cores)

# extract cellbender cleaned counts and perform initial cell-level filtering
process_seurat_object <- function(sample_type) {
  sct_path <- paste0(sample_type, "_sct.rds")
  
  # check if the SCT object already exists
  if (file.exists(sct_path)) {
    message("Preprocessed SCT object already exists.")
    return(invisible(NULL))
  }

  # load the merged object
  file_path <- paste0(sample_type, "_dual_seurat.rds")
  if (!file.exists(file_path)) {
    stop(paste("The merged Seurat object", file_path, "does not exist."))
  }
  merged_seurat_obj <- readRDS(file_path)
  
  # extract the counts 
  counts <- GetAssayData(object = merged_seurat_obj, assay = "RNA", layer = "counts")
    
  # create a new object 
  clean_seurat_obj <- CreateSeuratObject(counts = counts, meta.data = merged_seurat_obj@meta.data)
  
  # add metadata indicating the sample type
  clean_seurat_obj$sample_type <- sample_type
  
  # calculate percentage of mitochondrial genes
  clean_seurat_obj$percent.mito <- PercentageFeatureSet(clean_seurat_obj, pattern = "^mt-")
  
  # filter ribosomal protein genes
  rpS_gene <- grep(pattern = "^Rps", x = rownames(clean_seurat_obj), value = TRUE)
  rpL_gene <- grep(pattern = "^Rpl", x = rownames(clean_seurat_obj), value = TRUE)
  ps_gene <- grep(pattern = "-ps", x = rownames(clean_seurat_obj), value = TRUE)
  rpS_gene <- setdiff(rpS_gene, ps_gene)
  rpL_gene <- setdiff(rpL_gene, ps_gene)
  
  # calculate percentage of ribosomal protein genes
  clean_seurat_obj$percent.rpS <- PercentageFeatureSet(clean_seurat_obj, features = rpS_gene)
  clean_seurat_obj$percent.rpL <- PercentageFeatureSet(clean_seurat_obj, features = rpL_gene)
  clean_seurat_obj$percent.ps <- PercentageFeatureSet(clean_seurat_obj, features = ps_gene)
  
  # perform initial cell-level filtering
  clean_seurat_obj <- subset(clean_seurat_obj, subset = percent.mito <= 5 & nCount_RNA >= 500 & nCount_RNA <= 20000)
  clean_seurat_obj <- subset(clean_seurat_obj, subset = percent.rpL <= 5 & percent.rpS <= 5 & percent.ps <= 5)

  # save the cleaned object
  saveRDS(clean_seurat_obj, file = paste0(sample_type, "_cell_filtered.rds"))
  
  # split the object by samples
  clean_seurat_obj[["RNA"]] <- split(clean_seurat_obj[["RNA"]], f = clean_seurat_obj@meta.data$orig.ident)

  # Preprocessing steps
  # SCTransform normalization
  seurat_obj_sct <- SCTransform(
    clean_seurat_obj,
    vars.to.regress = c("nCount_RNA", "percent.mito", "percent.rpL", "percent.rpS", "percent.ps")
  )
  seurat_obj_sct <- RunPCA(seurat_obj_sct, reduction.name = "pca_sct", reduction.key = "PCASCT_")
  saveRDS(seurat_obj_sct, file = paste0(sample_type, "_sct.rds"))

  # clean up 
  rm(clean_seurat_obj, merged_seurat_obj)  
  gc()

  return(seurat_obj_sct)
}

# perform integration
perform_integration <- function(sample_type) {
  sct_path <- paste0(sample_type, "_sct.rds")
  integrated_path <- paste0(sample_type, "_sct_harmony.rds")
  
  if (!file.exists(sct_path)) {
    message("Preprocessed SCT object not found. Running process_seurat_object.")
    process_seurat_object(sample_type)
  }
  
  # load the preprocessed object 
  seurat_obj_sct <- readRDS(sct_path)
  
  # Harmony integration with SCT
  if (!"harmony_sct" %in% names(seurat_obj_sct@reductions)) {
    DefaultAssay(seurat_obj_sct) <- 'SCT'
    seurat_obj_sct <- IntegrateLayers(
      object = seurat_obj_sct,
      method = HarmonyIntegration,
      normalization.method = "SCT",
      orig.reduction = "pca_sct",
      new.reduction = "harmony_sct"
    )
  }
   
  # save the integrated data 
  saveRDS(seurat_obj_sct, file = integrated_path)

  return(seurat_obj_sct)
}

# find clusters and plot UMAP
run_clustering_and_plot <- function(seurat_obj, reduction_name, dims = 1:30, resolution = 0.5, algorithm = 1) {
  # ensure UMAP and clustering with integrated or unintegrated data
  seurat_obj <- FindNeighbors(seurat_obj, reduction = reduction_name, dims = dims)
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution, algorithm = algorithm)
  seurat_obj <- RunUMAP(seurat_obj, reduction = reduction_name, dims = dims, reduction.name = paste0("umap_", reduction_name))
  
  p1 <- DimPlot(seurat_obj, reduction = paste0("umap_", reduction_name), group.by = c("orig.ident","seurat_clusters") + ggtitle(paste("UMAP plot using", reduction_name, "reduction - Grouped by Clusters"))
  
  return(list(seurat_obj = seurat_obj, plot_clusters = p1))
}

# plotting and metrics
perform_plotting_and_metrics <- function(sample_type, output_dir) {
  sct_path <- paste0(sample_type, "_sct.rds")
  integrated_path <- paste0(sample_type, "_sct_harmony.rds")
    
  seurat_obj_integrated <- readRDS(integrated_path)
  
  pdf(file = file.path(output_dir, paste0(sample_type, "_dimplots.pdf")), width = 20, height = 10)

  # UMAP and clustering for unintegrated (PCA)
  pca_clustering_results <- run_clustering_and_plot(seurat_obj_sct, "pca_sct", dims = 1:30, resolution = 0.5, algorithm = 1)
  seurat_obj_sct <- pca_clustering_results$seurat_obj
  print(pca_clustering_results$plot_clusters)

  # UMAP and clustering for Harmony integrated data
  integration_results <- run_clustering_and_plot(seurat_obj_integrated, "harmony_sct", dims = 1:30, resolution = 0.5, algorithm = 1)
  seurat_obj_integrated <- integration_results$seurat_obj
  print(integration_results$plot_clusters)
  
  dev.off()

  # save the updated objects
  saveRDS(seurat_obj_integrated, file = integrated_path)
}

# run analysis
args <- commandArgs(trailingOnly = TRUE)
sample_type <- args[1]
output_dir <- "/path/to/folder/"  
perform_integration(sample_type)
perform_plotting_and_metrics(sample_type, output_dir)
