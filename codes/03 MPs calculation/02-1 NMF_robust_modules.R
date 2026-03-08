#!/usr/bin/env Rscript
###############################################################################
##  02_NMF_robust_modules.R
##
##  Perform non-negative matrix factorisation (NMF) on each GSM dataset
##  to identify robust gene modules per major cell type.
##
##  Author: wzy
##  Date  : 2025
###############################################################################

## ------------------------- PACKAGES & GLOBALS ------------------------------ ##
suppressPackageStartupMessages({
  library(NMF)
  library(dplyr)
  library(Seurat)
})

## NMF parameters
nrand       <- 3     # (kept for future randomisation tests)
nbin        <- 25    # (kept for future binning tests)
range_rank  <- 3:15  # candidate ranks to test
gmin        <- 5     # minimal gene number per module
pval_thresh <- 0.1   # (kept for downstream significance filtering)
seed        <- 777   # reproducibility

## ------------------------- HELPER FUNCTIONS -------------------------------- ##

##' Extract gene modules from an NMF fit
##'
##' @param res   NMF object (e.g. returned by nmf::nmf())
##' @param gmin  integer, minimal genes required per module
##' @return      named list of character vectors (gene names), or NA if failed
NMFToModules <- function(res, gmin = 5) {
  scores <- basis(res)      # genes × factors
  coefs  <- coefficients(res)  # factors × cells
  
  ## helper: identify top genes for each factor
  getTopGenes <- function(mat, gmin) {
    nr <- nrow(mat)
    nc <- ncol(mat)
    ranks_x <- apply(mat, 2, function(x) rank(-x, ties.method = "min"))
    ranks_y <- apply(mat, 2, rank)
    
    modules <- vector("list", nc)
    for (i in seq_len(nc)) {
      top <- which(ranks_x[, i] == 1)
      modules[[i]] <- names(top)[top]
    }
    len <- lengths(modules)
    modules <- modules[len >= gmin]
    if (length(modules) < 1) return(NULL)
    names(modules) <- paste0("m_", modules[[1]][1])
    modules
  }
  
  modules <- getTopGenes(scores, gmin)
  if (is.null(modules)) return(NA)
  modules
}

##' Basic Seurat preprocessing for a given cell type
##'
##' @param obj      Seurat object (full dataset)
##' @param celltype character, cell-type name to subset
##' @return         Seurat object (subset & processed)
preprocessData <- function(obj, celltype) {
  obj <- NormalizeData(obj) %>%
    FindVariableFeatures(nfeatures = 5000) %>%
    ScaleData(do.center = TRUE)
  subset(obj, subset = cell_type == celltype)
}

##' Run NMF across candidate ranks and pick the first rank that yields
##' exactly the same number of modules as the rank itself (heuristic stability).
##'
##' @param expr.list Seurat object (already subsetted to one cell type)
##' @return          named list of gene modules, or NA if insufficient modules
NMFdata <- function(expr.list) {
  mat <- as.matrix(GetAssayData(expr.list, slot = "scale.data"))
  mat <- mat[VariableFeatures(expr.list), ]
  mat[mat < 0] <- 0
  mat <- mat[rowSums(mat) > 0, ]
  
  modules.list <- NA
  for (r in rev(range_rank)) {          # try high → low
    set.seed(seed)
    fit <- nmf(mat, rank = r, nrun = 10, method = "nsNMF")
    tmp <- NMFToModules(fit, gmin = gmin)
    cat(sprintf("Rank %2d -> %2d modules\n", r, length(tmp)))
    if (identical(length(tmp), r)) {
      modules.list <- tmp
      break
    }
  }
  if (length(modules.list) >= 3) modules.list else NA
}

##' Identify cell types with >50 cells
##'
##' @param res.file Seurat object
##' @return         character vector of valid cell-type names
cellFind <- function(res.file) {
  meta <- res.file@meta.data
  tab  <- table(meta$cell_type)
  names(tab)[tab > 50]
}

## ------------------------- MAIN PIPELINE ----------------------------------- ##

data_dir  <- "/bioccXYJ/Teams_ZhangCL/pyy/pancreatic_cancer_scRNA/MP/prepare_data/"
out_dir   <- "."                       # change if you want separate folder
if (!dir.exists(out_dir)) dir.create(out_dir)

file_list <- list.files(data_dir, pattern = "^GSM.*\\.rds$")
message(sprintf("Found %d GSM files", length(file_list)))

for (f in file_list) {
  message("Processing: ", f)
  res.file <- readRDS(file.path(data_dir, f))
  
  cell_types <- cellFind(res.file)
  gsm_module <- list()
  
  for (ct in cell_types) {
    obj  <- preprocessData(res.file, celltype = ct)
    mods <- NMFdata(obj)
    if (!is.na(mods[[1]])) {
      names(mods) <- paste0("m_", seq_along(mods))  # ensure unique names
      gsm_module[[ct]] <- mods
    }
  }
  
  out_file <- file.path(out_dir, paste0("module_", f))
  saveRDS(gsm_module, file = out_file)
  message("Saved: ", out_file)
}

message("All done.")