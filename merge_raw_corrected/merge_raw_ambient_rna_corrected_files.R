"""
Author: Saranya Canchi
Date Created: 2024-08
Description:
- Merges CellBender and Cell Ranger outputs, constructs a dual-assay Seurat object, and saves the integrated dataset.
- Logs key memory usage and system information at each processing step to support resource tracking and reproducibility.
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
library(pryr)
library(parallel)

# get num_cores
num_cores <- parallel::detectCores()
message("Number of available cores: ", num_cores)

# get the dataset prefix 
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 1) {
  stop("Please provide a dataset prefix.")
}
dataset_prefix <- args[1]

# log memory usage
log_memory <- function(step, dataset_prefix, start_time) {
  mem_used <- pryr::mem_used()
  elapsed_time <- Sys.time() - start_time
  cat(paste("Step:", step, "for", dataset_prefix, "\n",
            "Memory used:", format(mem_used, units = "auto"), "\n",
            "Elapsed time:", elapsed_time, "seconds\n",
            "Timestamp:", Sys.time(), "\n\n"),
      file = paste0("slurm_logs/", dataset_prefix, "_memory_usage.log"), append = TRUE)
}

# process dataset and log
process_dataset <- function(dataset_prefix) {
  dual_seurat_output_file <- paste0(dataset_prefix, "_dual_seurat.rds")

  # start timer
  start_time <- Sys.time()
  
  # check if merged output exists
  if (!file.exists(dual_seurat_output_file)) {
    message("Merged output does not exist. Proceeding for dataset: ", dataset_prefix)

    # log initial memory usage
    log_memory("Initial state", dataset_prefix, start_time)
    
    # load or create merged cellbender output
    cellbender_output_file <- paste0(dataset_prefix, "_cellbender_merged.rds")
    if (file.exists(cellbender_output_file)) {
      message("Loading existing CellBender output for dataset: ", dataset_prefix)
      cellbender_merged <- readRDS(cellbender_output_file)
    } else {
      message("Creating CellBender output for dataset: ", dataset_prefix)
      cellbender_merged <- Read_CellBender_h5_Multi_Directory(
        base_path = paste0("/path/to/folder/", dataset_prefix),
        custom_name = "_cellbender.h5",
        merge = TRUE,
        parallel = TRUE,
        num_cores = num_cores
      )
      saveRDS(cellbender_merged, file = cellbender_output_file)
      log_memory("After merge cellbender outs", dataset_prefix, start_time)
    }

    # load or create merged cellranger output
    cellranger_output_file <- paste0(dataset_prefix, "_cellranger_merged.rds")
    if (file.exists(cellranger_output_file)) {
      message("Loading existing CellRanger output for dataset: ", dataset_prefix)
      cellranger_merged <- readRDS(cellranger_output_file)
    } else {
      message("Creating CellRanger output for dataset: ", dataset_prefix)
      cellranger_merged <- Read10X_h5_Multi_Directory(
        base_path = paste0("/path/to/folder/", dataset_prefix),
        h5_filename = "filtered_feature_bc_matrix.h5",
        merge = TRUE,
        parallel = TRUE,
        num_cores = num_cores
      )
      saveRDS(cellranger_merged, file = cellranger_output_file)
      log_memory("After merge cellranger outs", dataset_prefix, start_time)
    }

    # create dual assay seurat object
    message("Creating dual Seurat object for dataset: ", dataset_prefix)
    dual_seurat <- Create_CellBender_Merged_Seurat(
      raw_cell_bender_matrix = cellbender_merged,
      raw_counts_matrix = cellranger_merged,
      raw_assay_name = "RAW",
      min_cells = 5,
      min_features = 200
    )
    saveRDS(dual_seurat, file = dual_seurat_output_file)
    log_memory("After create dual assay seurat object", dataset_prefix, start_time)

    message("Finished processing for dataset: ", dataset_prefix)
  } else {
    message("Skipping dataset as output already exists: ", dataset_prefix)
  }
}

# log system information 
sys_info <- Sys.info()
cat(paste("System Information:\n",
          "Sysname:", sys_info['sysname'], "\n",
          "Release:", sys_info['release'], "\n",
          "Version:", sys_info['version'], "\n",
          "Machine:", sys_info['machine'], "\n\n"),
    file = paste0("slurm_logs/", dataset_prefix, "_system_info.log"), append = TRUE)

# run the processing function 
tryCatch({
  process_dataset(dataset_prefix)
}, error = function(e) {
  cat("Error processing dataset", dataset_prefix, ":", e$message, "\n")
})

# clean up
gc()
