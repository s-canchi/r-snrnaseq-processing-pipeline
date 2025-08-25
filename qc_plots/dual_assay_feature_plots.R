"""
Author: Saranya Canchi
Date Created: 2024-08
Description:
- Performs normalization, dimensionality reduction, clustering, and UMAP for both RNA and RAW assays in dual-assay Seurat object.
- Generates and saves feature plots for neuroscience-relevant marker genes.
"""

# load libraries
library(ggplot2)
library(dplyr)
library(magrittr)
library(patchwork)
library(viridis)
library(Seurat)
library(scCustomize)
library(qs)

# log memory usage
print_memory_usage <- function(stage) {
  memory <- pryr::mem_used()
  cat(stage, ": Memory usage = ", format(memory, units = "auto"), "\n")
}

# perform analysis steps (Normalization to UMAP) 
process_assay <- function(seurat_obj, assay_name) {
  # set active assay
  DefaultAssay(seurat_obj) <- assay_name

  # normalize data
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  
  # find variable features
  seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
  
  # scale data
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  
  # run PCA
  seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
  
  # find neighbors
  seurat_obj <- FindNeighbors(seurat_obj, verbose = FALSE)
  
  # find clusters
  seurat_obj <- FindClusters(seurat_obj, verbose = FALSE)
  
  # run UMAP
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:10, verbose = FALSE)

  print_memory_usage(paste("After processing", assay_name, "assay"))

  return(seurat_obj)
}

# generate feature plots
generate_feature_plots <- function(seurat_obj, markers, output_pdf) {
  pdf(output_pdf, width = 10, height = 8)
  
  # loop through each marker 
  for (marker_info in markers) {
    cell_type <- marker_info$cell_type
    marker <- marker_info$gene
    title <- paste(cell_type, ':', marker)
    
    print(paste("Plotting feature:", title))
    
    # generate the dual assay feature plot
    feature_plot <- FeaturePlot_DualAssay(seurat_obj, features = marker, assay1 = "RAW", assay2 = "RNA")
    
    # add the title using patchwork to combine with an empty plot with the title
    title_plot <- ggplot() + 
                  ggtitle(title) + 
                  theme_void() +
                  theme(plot.title = element_text(hjust = 0.5))  
    
    # combine the title plot and the feature plot using patchwork
    combined_plot <- title_plot / feature_plot + plot_layout(heights = c(1, 10))
    
    print(combined_plot)
  }
  
  dev.off()
}

# define a list of markers with cell types
marker_groups <- list(
  Neurons = c("Rbfox3", "Map2", "Syt1", "Snap25"),
  Astrocytes = c("Gfap", "Aldh1l1", "S100b"),
  Oligodendrocytes = c("Mbp", "Plp1", "Cnp"),
  Microglia = c("Aif1", "Tmem119", "Cx3cr1"),
  "Endothelial Cells" = c("Cldn5", "Pecam1", "Flt1"),
  "Exc. Neurons" = c("Slc17a7", "Slc17a6", "Camk2a"),
  "Inh. Neurons" = c("Gad1", "Gad2", "Pvalb", "Sst", "Vip")
)

mouse_brain_markers <- unlist(
  lapply(names(marker_groups), function(celltype)
    lapply(marker_groups[[celltype]], function(gene)
      list(cell_type = celltype, gene = gene)
    )
  ), recursive = FALSE
)

# main args
dataset_prefixes <- c("prefix1", "prefix2")
input_dir <- "/path/to/folder/"
output_dir <- "/path/to/folder/"

# set wd
setwd(input_dir)

# loop through each dataset prefix 
for (prefix in dataset_prefixes) {
  seurat_obj_path <- paste0(input_dir, "/", prefix, "_dual_seurat.rds")
  output_path <- paste0(input_dir, "/", prefix, "_dual_seurat_normalized.rds")
  feature_plots_pdf <- paste0(output_dir, "/", prefix, "_feature_plots.pdf")
  
  # check if output file already exists
  if (file.exists(output_path)) {
    cat("Output file", output_path, "already exists. Skipping processing for", prefix, "\n")
  } else {
    # load object
    seurat_obj <- readRDS(seurat_obj_path)
    print_memory_usage(paste("After loading Seurat object for", prefix))
    
    # process RNA assay
    seurat_obj <- process_assay(seurat_obj, "RNA")
    
    # process RAW assay
    seurat_obj <- process_assay(seurat_obj, "RAW")
    
    # save the processed object
    saveRDS(seurat_obj, output_path)
    print_memory_usage(paste("After saving processed object for", prefix))
    
    cat("Processing of", prefix, "assays completed successfully.\n")
  }
  
  # generate and save the feature plots
  generate_feature_plots(seurat_obj, mouse_brain_markers, feature_plots_pdf)
}
  
