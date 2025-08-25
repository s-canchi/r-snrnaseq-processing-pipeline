"""
Author: Saranya Canchi
Date Created: 2024-08
Description:
- Integrates doublet detection metadata into Seurat objects and flags high-confidence doublets for removal.
- Generates UMAP plots and summary statistics to visualize doublet identification and the effects of filtering.
"""

# load libraries
library(Seurat)
library(dplyr)
library(ggplot2)

# set paths
setwd("/path/to/folder/")
output_dir <- "/path/to/folder/"
intermediate_dir <- file.path(getwd(), "doublet_intermediate_files")
params_file <- file.path(intermediate_dir, "doubletfinder_parameters.csv")

# load metadata list
load_metadata_list <- function(sample_type, intermediate_dir) {
  metadata_file <- file.path(intermediate_dir, sprintf("%s_metadata_list.rds", sample_type))
  if (!file.exists(metadata_file)) {
    stop("Metadata list file not found.")
  }
  metadata_list <- readRDS(metadata_file)
  combined_metadata <- bind_rows(metadata_list)
  print(head(combined_metadata,10))

  # check the combined metadata has the correct columns
  required_cols <- c("doublet_probability", "DF.classification_initial", "DF.classification_adjusted")
  if (!all(required_cols %in% colnames(combined_metadata))) {
    stop("Combined metadata does not have the expected columns: ", paste(required_cols, collapse = ", "))
  }

  return(combined_metadata)
}

# define sample types
sample_types <- c("hpc", "pfc")

# perform doublet detection and QC plotting
run_doublet_qc <- function(sample_type, intermediate_dir, output_dir) {
  orig_filename <- sprintf("%s_sct_harmony.rds", sample_type)
  seurat_obj <- readRDS(orig_filename)
  
  # load combined metadata
  combined_metadata <- load_metadata_list(sample_type, intermediate_dir)
  
  # check if combined metadata is present
  if (nrow(combined_metadata) == 0) {
    stop("No metadata found. Please ensure that the intermediate files are correctly generated.")
  }
  
  # check cell barcodes of combined and original metadata 
  cat("First few rows of combined metadata:\n")
  print(head(combined_metadata))
  
  cat("First few rows of original Seurat object metadata cell barcodes:\n")
  print(head(rownames(seurat_obj@meta.data)))

  # ensure the cell barcodes are consistent
  combined_metadata <- combined_metadata %>%
    dplyr::mutate(cell_barcodes = rownames(.))
  seurat_obj@meta.data <- seurat_obj@meta.data %>%
    dplyr::mutate(cell_barcodes = rownames(.))

  # ensure overlap between cell barcodes
  intersect_barcodes <- intersect(seurat_obj@meta.data$cell_barcodes, combined_metadata$cell_barcodes)
  if (length(intersect_barcodes) == 0) {
    stop("No cell overlap between new meta data and Seurat object")
  }

  cat("Number of overlapping barcodes:", length(intersect_barcodes), "\n")

  # merge combined metadata with original Seurat object metadata
  combined_metadata <- combined_metadata %>%
    dplyr::select(cell_barcodes, doublet_probability, DF.classification_initial, DF.classification_adjusted)

  seurat_obj@meta.data <- seurat_obj@meta.data %>%
    dplyr::left_join(combined_metadata, by = "cell_barcodes")
  rownames(seurat_obj@meta.data) <- seurat_obj@meta.data$cell_barcodes
  
  # set the doublet classification and doublet probability columns
  seurat_obj$doublet_probability <- seurat_obj@meta.data$doublet_probability
  seurat_obj$Doublets <- seurat_obj@meta.data$DF.classification_adjusted

  print(head(seurat_obj@meta.data, 10))
 
  # high-resolution clustering
  seurat_obj <- FindClusters(seurat_obj, graph.name = "harmony_sct_snn", resolution = 1.5, cluster.name = "doublet_cleanup")
  
  # identify clusters with high doublet proportions
  cluster_doublet_proportions <- seurat_obj@meta.data %>%
    group_by(doublet_cleanup) %>%
    summarize(doublet_prop = mean(Doublets == "Doublet"))

  clusters_to_flag <- cluster_doublet_proportions %>%
    filter(doublet_prop > 0.70) %>%
    pull(doublet_cleanup)
  
  # flag nuclei as doublets if they belong to a high-doublet cluster
  seurat_obj$final_doublet_label <- seurat_obj$Doublets == "Doublet" | seurat_obj$doublet_cleanup %in% clusters_to_flag

  # filter out doublets
  filtered_seurat_obj <- subset(seurat_obj, final_doublet_label == FALSE)

  # generate qc plots 
  pdf_file <- file.path(output_dir, sprintf("%s_doublet_qc_plots.pdf", sample_type))
  pdf(file = pdf_file, width = 16, height = 8.5)
  
  # doublet confidence plot
  seurat_obj$doublet_confidence <- "Singlet"
  seurat_obj$doublet_confidence[grepl("Doublet", seurat_obj$Doublets) & seurat_obj$doublet_probability <= 0.5] <- "Doublet_LowConf"
  seurat_obj$doublet_confidence[grepl("Doublet", seurat_obj$Doublets) & seurat_obj$doublet_probability > 0.5] <- "Doublet_HighConf"
  
  print(DimPlot(seurat_obj, reduction = "umap_sct_harmony", group.by = "doublet_confidence", cols = c("red", "gray", "yellow")) + 
        ggtitle(sprintf("%s UMAP - Doublet Confidence", sample_type)))

  # pre-filtering UMAP plot
  print(DimPlot(seurat_obj, reduction = "umap_sct_harmony", group.by = c("doublet_cleanup","Doublets")) + 
  scale_color_manual(values = c("Singlet" = "gray", "Doublet" = "red")) +
  ggtitle(sprintf("%s UMAP - Before Filtering", sample_type)))
  
  # post-filtering UMAP plot
  print(DimPlot(filtered_seurat_obj, reduction = "umap_sct_harmony", group.by = c("doublet_cleanup","Doublets")) + 
  scale_color_manual(values = c("Singlet" = "gray", "Doublet" = "red")) +
  ggtitle(sprintf("%s UMAP - After Filtering", sample_type)))  

  dev.off()
  message(sprintf("All plots saved to %s", pdf_file))

  # save filtered object
  saveRDS(filtered_seurat_obj, sprintf("%s_doublet_filtered.rds", sample_type))

  # clean up
  gc()

  return(pdf_file)
}

# run analysis
for (sample_type in sample_types) {
  run_doublet_qc(sample_type, intermediate_dir, output_dir)
}