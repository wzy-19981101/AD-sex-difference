#!/usr/bin/env Rscript
###############################################################################
##  DEG calculation.R
##
##  Differential Expression Analysis for scRNA-seq data (AD vs normal,
##  female vs male) across multiple datasets.
##
##  Author: wzy
##  Date  : 2025
###############################################################################

## ------------------------- PACKAGES & GLOBALS ------------------------------ ##

library(Seurat)
library(SeuratObject)
library(readr)
library(harmony)

# ----------------------------------------------------------------------------
# Configuration: Set base directories
# ----------------------------------------------------------------------------
base_dir <- "/bioccXYJ/Teams_ZhangCL/wzy/ADscRNAseq"
metadata_path <- file.path(base_dir, "demarkers_sex_disease", "allmetadata.rds")

# ROSMAP data
rosmap_data_dir <- file.path(base_dir, "data_celltype", "ROSMAP427")
rosmap_output_dir <- file.path(base_dir, "demarkers_sex_disease", "new427result")

# GEO data
geo_base_dir <- "/bioccXYJ/Teams_ZhangCL/wzy/zyx_mp/output/GEO"
geo_output_base <- file.path(base_dir, "demarkers_sex_disease", "GEO")

# Create output directories if they don't exist
dir.create(rosmap_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(geo_output_base, recursive = TRUE, showWarnings = FALSE)

# Load subject metadata (contains sex, disease, age, region, apoe)
message("Loading metadata...")
msex <- readRDS(metadata_path)

# ----------------------------------------------------------------------------
# Helper function: Run differential expression for a given subset
# ----------------------------------------------------------------------------
#' Run differential expression analysis using Seurat's FindMarkers
#'
#' @param obj Seurat object subset (e.g., females only)
#' @param ident.1 First identity class (e.g., "AD")
#' @param ident.2 Second identity class (e.g., "normal")
#' @param group.by Metadata column to define identities (e.g., "disease")
#' @param latent.vars Character vector of latent variables to regress out
#' @param logfc.threshold Log fold change threshold for markers
#' @param output_path Full path to save the result RDS file
#' @param min_cells Minimum number of cells required in each group to proceed
#'
#' @return Invisibly returns the markers data frame, or NULL if skipped
run_de_analysis <- function(obj, ident.1, ident.2, group.by,
                            latent.vars = c("age", "region", "apoe_genotype"),
                            logfc.threshold = 0.25,
                            output_path,
                            min_cells = 200) {
  
  # Check if group.by column exists
  if (!(group.by %in% colnames(obj@meta.data))) {
    warning("Column '", group.by, "' not found in metadata. Skipping.")
    return(NULL)
  }
  
  # Count cells per group
  cell_counts <- table(obj@meta.data[[group.by]])
  if (length(cell_counts) < 2) {
    warning("Not enough groups in '", group.by, "'. Skipping.")
    return(NULL)
  }
  
  # Ensure both groups have at least min_cells
  counts_needed <- c(ident.1, ident.2)
  if (!all(counts_needed %in% names(cell_counts))) {
    warning("Required groups (", paste(counts_needed, collapse = ", "), 
            ") not both present. Skipping.")
    return(NULL)
  }
  if (any(cell_counts[counts_needed] < min_cells)) {
    message("Skipping due to insufficient cells: ", 
            paste(counts_needed, cell_counts[counts_needed], sep = "=", collapse = ", "))
    return(NULL)
  }
  
  # check apoe_genotype
  latent_vars_adj <- latent.vars
  if ("apoe_genotype" %in% latent_vars_adj) {
    
    if (!"apoe_genotype" %in% colnames(obj@meta.data)) {
      warning("Column 'apoe_genotype' not found in metadata. Removing from latent.vars.")
      latent_vars_adj <- setdiff(latent_vars_adj, "apoe_genotype")
    } else {
     
      if (all(is.na(obj@meta.data$apoe_genotype))) {
        message("apoe_genotype is all NA. Removing it from latent.vars.")
        latent_vars_adj <- setdiff(latent_vars_adj, "apoe_genotype")
      }
    }
  }
  # --------------------------------------------------------------------------
  
  # Set identities
  Idents(obj) <- obj@meta.data[[group.by]]
  
  # Run FindMarkers 
  markers <- FindMarkers(obj,
                         test.use = "wilcox",
                         ident.1 = ident.1,
                         ident.2 = ident.2,
                         latent.vars = latent_vars_adj,
                         pseudocount.use = 1,
                         logfc.threshold = logfc.threshold)
  
  # Save result
  saveRDS(markers, file = output_path)
  message("Saved: ", output_path)
  
  invisible(markers)
}

# ----------------------------------------------------------------------------
# Process ROSMAP datasets
# ----------------------------------------------------------------------------
message("\n=== Processing ROSMAP datasets ===")

# List cell type files 
rosmap_files <- list.files(rosmap_data_dir, pattern = "\\.rds$", full.names = TRUE)


for (file_path in rosmap_files) {
  celltype <- basename(file_path)
  message("\nProcessing ROSMAP cell type: ", celltype)
  
  # Read data
  data <- tryCatch(readRDS(file_path), error = function(e) {
    warning("Failed to read ", file_path, ": ", e$message)
    return(NULL)
  })
  if (is.null(data)) next
  
  # Add subject metadata (sex, disease, etc.)
  subject_df <- data.frame(subject = data@meta.data$subject, stringsAsFactors = FALSE)
  meta_to_add <- msex[match(subject_df$subject, msex$subject), ]
  data <- AddMetaData(data, metadata = meta_to_add)
  
  # Standard preprocessing
  data <- NormalizeData(data)
  data <- FindVariableFeatures(data)
  data <- ScaleData(data)
  data <- RunPCA(data)
  
  # Subset by sex and disease
  data_F <- subset(data, subset = msex == "Female")
  data_M <- subset(data, subset = msex == "Male")
  data_AD <- subset(data, subset = disease == "AD")
  
  # Define output base name
  out_base <- file.path(rosmap_output_dir, paste0(celltype, "_markers"))
  
  # Female: AD vs normal
  run_de_analysis(data_F,
                  ident.1 = "AD", ident.2 = "normal",
                  group.by = "disease",
                  output_path = paste0(out_base, "_F.rds"))
  
  # Male: AD vs normal
  run_de_analysis(data_M,
                  ident.1 = "AD", ident.2 = "normal",
                  group.by = "disease",
                  output_path = paste0(out_base, "_M.rds"))
  
  # AD: Female vs Male
  run_de_analysis(data_AD,
                  ident.1 = "Female", ident.2 = "Male",
                  group.by = "msex",
                  output_path = paste0(out_base, "_AD.rds"))
}

# ----------------------------------------------------------------------------
# Process GEO datasets
# ----------------------------------------------------------------------------
message("\n=== Processing GEO datasets ===")

geo_datasets <- list.files(geo_base_dir)
# Select only the 8th dataset? Original code used datasets[8]. Adjust as needed.
# Here we loop over all datasets found; you can filter if necessary.
# geo_datasets <- geo_datasets[8]  # Uncomment to process only the 8th

for (dataset in geo_datasets) {
  message("\nProcessing GEO dataset: ", dataset)
  
  # Create output subdirectory
  out_dir <- file.path(geo_output_base, dataset)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # List RDS files in the dataset folder
  dataset_path <- file.path(geo_base_dir, dataset)
  rds_files <- list.files(dataset_path, pattern = "\\.rds$", full.names = TRUE)
  
  for (file_path in rds_files) {
    celltype <- basename(file_path)
    message("  Processing cell type: ", celltype)
    
    # Read data
    data <- tryCatch(readRDS(file_path), error = function(e) {
      warning("Failed to read ", file_path, ": ", e$message)
      return(NULL)
    })
    if (is.null(data)) next
    
    # Add subject metadata
    subject_df <- data.frame(subject = data@meta.data$orig.ident, stringsAsFactors = FALSE)
    meta_to_add <- msex[match(subject_df$subject, msex$subject), ]
    data <- AddMetaData(data, metadata = meta_to_add)
    
    # Preprocessing
    data <- NormalizeData(data)
    data <- FindVariableFeatures(data)
    data <- ScaleData(data)
    data <- RunPCA(data)
    
    # Determine disease column name (either "disease" or "cancer_type")
    disease_col <- if ("disease" %in% colnames(data@meta.data)) "disease" else "cancer_type"
    
    # Subsets
    data_F <- subset(data, subset = msex == "Female")
    data_M <- subset(data, subset = msex == "Male")
    data_AD <- subset(data, subset = data@meta.data[[disease_col]] == "AD")
    
    # Output base
    out_base <- file.path(out_dir, paste0(celltype, "_markers"))
    
    # Female: AD vs normal
    run_de_analysis(data_F,
                    ident.1 = "AD", ident.2 = "normal",
                    group.by = disease_col,
                    logfc.threshold = 1,  # GEO uses stricter threshold
                    output_path = paste0(out_base, "_F.rds"))
    
    # Male: AD vs normal
    run_de_analysis(data_M,
                    ident.1 = "AD", ident.2 = "normal",
                    group.by = disease_col,
                    logfc.threshold = 1,
                    output_path = paste0(out_base, "_M.rds"))
    
    # AD: Female vs Male
    run_de_analysis(data_AD,
                    ident.1 = "Female", ident.2 = "Male",
                    group.by = "msex",
                    logfc.threshold = 1,
                    output_path = paste0(out_base, "_AD.rds"))
  }
}

message("\nAll processing completed.")