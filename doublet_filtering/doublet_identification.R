"""
Author: Saranya Canchi
Date Created: 2024-08
Description:
- Detects and annotates doublets in snRNA-seq samples using DoubletFinder.
- Saves results and doublet classification metadata for each sample.
"""

# load libraries
library(Seurat)
library(DoubletFinder)
library(ROCR)
library(dplyr)
library(ggplot2)
library(fs)

# set paths
setwd("/path/to/folder/")
intermediate_dir <- file.path(getwd(), "doublet_intermediate_files")
params_file <- file.path(intermediate_dir, "doubletfinder_parameters.csv")

# define filename pattern
filename_pattern <- function(sample_type) {
  sprintf("%s_sct_harmony.rds", sample_type)
}

# set doublet rate
doublet_rate <- 0.05  

# save DoubletFinder parameters
save_params <- function(sample_type, sample_id, optimal_pK, homotypic_proportion, expected_doublelets, adjusted_expected_doublelets) {
  params_data <- data.frame(
    sample_type = sample_type,
    sample_id = sample_id,
    optimal_pK = optimal_pK,
    homotypic_proportion = homotypic_proportion,
    expected_doublelets = expected_doublelets,
    adjusted_expected_doublelets = adjusted_expected_doublelets
  )
  
  # write the params to file
  if (!file.exists(params_file)) {
    write.csv(params_data, params_file, row.names = FALSE)
  } else {
    write.csv(params_data, params_file, row.names = FALSE, append = TRUE, col.names = FALSE)
  }
}

# load Seurat objects, identify doublets, filter, and generate QC plots
run_doubletfinder <- function(sample_type, intermediate_dir) {
  filename <- filename_pattern(sample_type)
  seurat_obj <- readRDS(filename)
  print(intermediate_dir)

  # get samples ids
  sample_ids <- as.vector(unique(seurat_obj$orig.ident))
  metadata_list <- list()

  # identify doublets for each sample
  for (sample_id in sample_ids) {
    print(paste("Processing:", sample_id))

    # check if the sample file exists 
    rds_file <- file.path(intermediate_dir, sprintf("%s_%s_doubletFinder.rds", sample_type, sample_id))
    if (file.exists(rds_file)) {
      message(sprintf("Sample %s already processed. Skipping...", sample_id))
      next
    }

    sample_subset <- subset(seurat_obj, orig.ident == sample_id)
    if (nrow(sample_subset@meta.data) == 0) {
      next  
    }

    # homotypic doublet proportion estimate based on specified clustering results
    homotypic_proportion <- modelHomotypic(Idents(sample_subset))  
    expected_doublelets <- round(doublet_rate * ncol(sample_subset))  
    adjusted_expected_doublelets <- round(expected_doublelets * (1 - homotypic_proportion))

    # initial pK identification (no ground-truth)
    sweep_results <- paramSweep(sample_subset, PCs = 1:30, sct = TRUE)
    sweep_stats <- summarizeSweep(sweep_results, GT = FALSE)
    bcmvn <- find.pK(sweep_stats)
    optimal_pK <- as.numeric(levels(bcmvn$pK)[bcmvn$pK][which.max(bcmvn$BCmetric)])

    print(sample_id)
    print(optimal_pK)
    print(homotypic_proportion)
    print(expected_doublelets)
    print(adjusted_expected_doublelets)

    # save parameters 
    save_params(sample_type, sample_id, optimal_pK, homotypic_proportion, expected_doublelets, adjusted_expected_doublelets)

    # initial DoubletFinder run
    sample_subset <- doubletFinder(sample_subset, PCs = 1:30, pN = 0.5, pK = optimal_pK, nExp = expected_doublelets, reuse.pANN = FALSE, sct = TRUE)
    
    # adjusted DoubletFinder run
    sample_subset <- doubletFinder(sample_subset, PCs = 1:30, pN = 0.5, pK = optimal_pK, nExp = adjusted_expected_doublelets, reuse.pANN = sprintf("pANN_0.5_%s_%s", optimal_pK, expected_doublelets), sct = TRUE)

    # rename columnnames
    meta_cols <- colnames(sample_subset@meta.data)
    pANN_col <- grep("^pANN_", meta_cols, value = TRUE)
    classification_cols <- grep("^DF.classifications_", meta_cols, value = TRUE)
    initial_classification_col <- classification_cols[order(match(classification_cols, meta_cols))[1]]
    adjusted_classification_col <- classification_cols[order(match(classification_cols, meta_cols))[2]]

   if (length(pANN_col) == 1 && length(initial_classification_col) == 1 && length(adjusted_classification_col) == 1) {
    colnames(sample_subset@meta.data)[meta_cols == pANN_col] <- "doublet_probability"
    colnames(sample_subset@meta.data)[meta_cols == initial_classification_col] <- "DF.classification_initial"
    colnames(sample_subset@meta.data)[meta_cols == adjusted_classification_col] <- "DF.classification_adjusted"
   } else {
    stop("Error: Expected metadata columns not found.")
   }

  # store metadata 
  metadata_list[[sample_id]] <- sample_subset@meta.data

  # save intermediate object
  saveRDS(sample_subset, rds_file)

  # clean up
  gc()

 }
  # save the metadata 
  saveRDS(metadata_list, file.path(intermediate_dir, sprintf("%s_metadata_list.rds", sample_type)))
}
  
# process all samples
args <- commandArgs(trailingOnly = TRUE)
sample_type <- args[1]
run_doubletfinder(sample_type, intermediate_dir)