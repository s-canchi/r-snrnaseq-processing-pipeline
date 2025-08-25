# snRNAseq Analysis 

This repository contains example R scripts for single-nucleus RNA-seq (snRNA-seq) data processing, quality control, doublet filtering, integration, and visualization, developed for neuroscience research. The scripts use Seurat and related packages and are organized by workflow stage.

## Repository Structure

- **cellranger/**  
  Scripts to generate count matrixes via CellRanger.

- **merge_raw_corrected/**  
  Combines Cell Ranger (raw) and CellBender (corrected) data to compute QC metrics and comparative plots.

- **qc_plots/**  
  Generates summary quality control visualizations and statistics for each sample/dataset.

- **doublet_filtering/**  
  Detects and removes doublets; produces visualizations of doublet status before and after filtering.

- **integration/**  
  Performs normalization, batch correction, and data integration.

## Notes

- Scripts are intended as modular workflow examples, not a complete automated pipeline.
- Paths and resource parameters are example values; adjust as needed for your local environment and data.



