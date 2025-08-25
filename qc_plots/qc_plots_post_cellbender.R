"""
Author: Saranya Canchi
Date Created: 2024-08
Description:
- Summarizes and compares key QC and count metrics before and after ambient RNA correction.
- Generates and saves summary plots and tables to visualize median metrics and gene-level changes across samples in each dataset.
"""

# load libraries
library(ggplot2)
library(dplyr)
library(reshape2)
library(tidyr)
library(grid)
library(magrittr)
library(patchwork)
library(viridis)
library(Seurat)
library(scCustomize)
library(qs)

# set wd
setwd("/path/to/folder/")

# list of datasets
datasets <- list(
  "prefix1" = "/path/to/folder/",
  "prefix2" = "/path/to/folder/"
)

# process each dataset
process_dataset <- function(prefix, dataset_path) {
  # define output file paths
  pdf_file <- paste0(prefix, "_median_qc_metrics_combined_barplot.pdf")
  csv_file <- paste0(prefix, "_feature_diff.csv")

  # check if output files already exist
  if (file.exists(pdf_file) && file.exists(csv_file)) {
    message(paste("Skipping dataset", prefix, "- output files already exist."))
    return()
  } else {
    message(paste("Processing dataset", prefix, "- output files do not exist."))
  }

  # load dataset
  dual_assay <- readRDS(dataset_path)
  message("Loaded dataset: ", dataset_path)

  # calculate diff pre/post CellBender
  dual_assay <- Add_CellBender_Diff(dual_assay, raw_assay_name = "RAW", cell_bender_assay_name = "RNA")

  # Part I: calculate and prep df for plotting median values of count, features, and diff values
  # calculate median stats
  median_stats <- Median_Stats(dual_assay, group_by_var="orig.ident", median_var=c("nCount_Diff","nFeature_Diff","nCount_RAW","nFeature_RAW"))

  # create a long format data frame
  medians_long <- melt(median_stats, id.vars = "orig.ident", variable.name = "Metric", value.name = "Median")

  # separate the metric names and types
  medians_long <- medians_long %>%
    mutate(Metric_Type = case_when(
      grepl("Count_RAW$", Metric) ~ "Raw",
      grepl("Feature_RAW$", Metric) ~ "Raw",
      grepl("Count_RNA$", Metric) ~ "Filtered",
      grepl("Feature_RNA$", Metric) ~ "Filtered",
      grepl("Diff$", Metric) ~ "Diff",  # Handle Diff metrics
      TRUE ~ NA_character_  # Default to NA for any other cases
    )) %>%
    mutate(Metric_Name = sub("Median_", "", Metric),
           Metric_Name = sub("_RAW", "", Metric_Name),
           Metric_Name = sub("_RNA", "", Metric_Name),
           Metric_Name = sub("_Diff", "", Metric_Name))

  # separate data for different plots accurately
  count_data <- filter(medians_long, Metric_Name == "nCount" & Metric_Type != "Diff")
  feature_data <- filter(medians_long, Metric_Name == "nFeature" & Metric_Type != "Diff")
  diff_data <- filter(medians_long, Metric_Type == "Diff")

  # calculate across-sample average values
  count_avg <- mean(count_data$Median, na.rm = TRUE)
  feature_avg <- mean(feature_data$Median, na.rm = TRUE)
  count_diff_avg <- mean(diff_data$Median[diff_data$Metric == "Median_nCount_Diff"], na.rm = TRUE)
  feature_diff_avg <- mean(diff_data$Median[diff_data$Metric == "Median_nFeature_Diff"], na.rm = TRUE)

  # Part II: calculate feature set difference
  # generate feature diff df
  feature_diff <- CellBender_Feature_Diff(dual_assay, raw_assay = "RAW", cell_bender_assay = "RNA")
  write.csv(feature_diff, paste0(prefix, "_feature_diff.csv"))

  # Part III: calculate QC differences between pre/post CellBender
  # calculate QC metrics pre CellBender
  raw_assay_obj <- Add_Cell_QC_Metrics(dual_assay, species="mouse", add_cell_cycle=FALSE, assay="RAW")

  # extract raw metadata
  raw_counts <- raw_assay_obj@meta.data
  raw_counts$type <- "Raw"

  # calculate QC metrics post CellBender
  filter_assay_obj <- Add_Cell_QC_Metrics(dual_assay, species="mouse", add_cell_cycle=FALSE, assay="RNA")

  # extract filter metadata
  filtered_counts <- filter_assay_obj@meta.data
  filtered_counts$type <- "Filtered"

  # combine the raw and filtered data frames
  combined_counts <- bind_rows(raw_counts, filtered_counts)

  # find column positions 
  start_col <- which(names(combined_counts) == "percent_mito")
  end_col <- which(names(combined_counts) == "percent_ieg")
  final_cols <- names(combined_counts)[start_col:end_col]

  # calculate medians for each sample and type
  medians <- combined_counts %>%
    group_by(orig.ident, type) %>%
    summarise(across(all_of(final_cols), \(x) median(x, na.rm = TRUE)), .groups = 'drop')

  # reshape the data frame to long format
  median_long <- melt(medians, id.vars = c("orig.ident", "type"), variable.name = "Metric", value.name = "Median")

  # create and save all plots to a single PDF file
  pdf(file = paste0(prefix, "_median_qc_metrics_cellbender.pdf"), width = 10, height = 8)

  # plot counts
  p1 <- ggplot(count_data, aes(x = orig.ident, y = Median, fill = Metric_Type)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.5), alpha = 0.7) +
    geom_hline(yintercept = count_avg, linetype = "dashed", color = "red") +
    annotate("text", x = -Inf, y = count_avg, label = round(count_avg, 2), vjust = -0.5, hjust = -0.1, color = "red") +
    labs(title = paste0("Median nCount for Each Sample (", prefix, " pre/post CellBender)"), x = "Sample", y = "Median Value") +
    scale_fill_manual(values = c("Raw" = "#1f77b4", "Filtered" = "#ff7f0e")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p1)

  # plot features
  p2 <- ggplot(feature_data, aes(x = orig.ident, y = Median, fill = Metric_Type)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.5), alpha = 0.7) +
    geom_hline(yintercept = feature_avg, linetype = "dashed", color = "red") +
    annotate("text", x = -Inf, y = feature_avg, label = round(feature_avg, 2), vjust = -0.5, hjust = -0.1, color = "red") +
    labs(title = paste0("Median nFeature for Each Sample (", prefix, " pre/post CellBender)"), x = "Sample", y = "Median Value") +
    scale_fill_manual(values = c("Raw" = "#1f77b4", "Filtered" = "#ff7f0e")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p2)

  # plot count diff
  p3 <- ggplot(filter(diff_data, Metric == "Median_nCount_Diff"), aes(x = orig.ident, y = Median)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.5), fill = "gray") +
    geom_hline(yintercept = count_diff_avg, linetype = "dashed", color = "red") +
    annotate("text", x = -Inf, y = count_diff_avg, label = round(count_diff_avg, 2), vjust = -0.5, hjust = -0.1, color = "red") +
    labs(title = paste0("Median nCount Diff for Each Sample (", prefix, " pre/post CellBender)"), x = "Sample", y = "Median Value") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p3)

  # plot feature diff
  p4 <- ggplot(filter(diff_data, Metric == "Median_nFeature_Diff"), aes(x = orig.ident, y = Median)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.5), fill = "gray") +
    geom_hline(yintercept = feature_diff_avg, linetype = "dashed", color = "red") +
    annotate("text", x = -Inf, y = feature_diff_avg, label = round(feature_diff_avg, 2), vjust = -0.5, hjust = -0.1, color = "red") +
    labs(title = paste0("Median nFeature Diff for Each Sample (", prefix, " pre/post CellBender)"), x = "Sample", y = "Median Value") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p4)

  # plot gene feature diff
  p5 <- CellBender_Diff_Plot(feature_diff, pct_diff_threshold=50, num_labels=30)
  print(p5)

  # loop through each unique metric and create a separate plot
  unique_metrics <- unique(median_long$Metric)
  for (metric in unique_metrics) {
    metric_data <- subset(median_long, Metric == metric)
    p <- ggplot(metric_data, aes(x = orig.ident, y = Median, fill = type)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.5), alpha = 0.7) +
      labs(title = paste("Median", metric, "for each sample pre/post CellBender"), x = "Sample", y = "Median Value") +
      scale_fill_manual(values = c("Raw" = "#1f77b4", "Filtered" = "#ff7f0e")) +  
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    print(p)
  }

  dev.off()
  
  # clean up
  gc()
}

# loop through the datasets 
for (prefix in names(datasets)) {
  process_dataset(prefix, datasets[[prefix]])
}

