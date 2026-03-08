#!/usr/bin/env Rscript
###############################################################################
##  02-3_filter_orthogonal_genes.R
##
##  From NMF modules and Leiden clusters:
##  1. Collect genes per cluster (top 50 % by frequency)
##  2. Remove shared genes between modules (keep higher frequency)
##  3. Save orthogonal gene lists (≥ 3 genes per module)
##
##  Author: wzy
##  Date  : 2025
##  Depends: base R
############################################################################--

## ------------------------- PARAMETERS ------------------------------------- ##
CELLTYPE      <- "oligodendrocytes"
ALL_MODULES   <- "D:/wuzhiyigongzuo/2024.9.12/nmf_zyy/program/all_modules.rds"
CLUSTER_FILE  <- paste0("D:/wuzhiyigongzuo/2024.9.12/nmf_zyy/program/NC/cluster_programs/", CELLTYPE, ".rds")
OUT_DIR       <- "D:/wuzhiyigongzuo/2024.9.12/nmf_zyy/program/NC/orthogonal_genes"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## ------------------------- LOAD DATA -------------------------------------- ##
all_modules  <- readRDS(ALL_MODULES)
module_genes <- all_modules[[CELLTYPE]]   # list: module_name -> genes
cluster_res  <- readRDS(CLUSTER_FILE)     # data.frame: cluster -> modules

## ------------------------- TOP 50 % BY FREQUENCY -------------------------- ##
MPs     <- list()   # final orthogonal genes
MPs_df  <- list()   # frequency table per module

for (cl in colnames(cluster_res)) {
  ## collect all genes in this cluster
  mods   <- na.omit(cluster_res[, cl])
  genes  <- unlist(module_genes[mods], use.names = FALSE)
  
  ## frequency table
  freq <- sort(table(genes), decreasing = TRUE)
  n_unique <- length(unique(freq))
  top50_th  <- unique(freq)[max(1, floor(n_unique * 0.5))]
  top_genes <- names(freq)[freq >= top50_th]
  
  MPs[[paste0("MP", cl)]] <- top_genes
  MPs_df[[paste0("MP", cl)]] <- data.frame(genes = names(freq),
                                           Freq  = as.numeric(freq),
                                           row.names = names(freq))
}

## ------------------------- REMOVE SHARED GENES ---------------------------- ##
module_names <- names(MPs_df)

for (i in seq_along(module_names)) {
  mod_i <- module_names[i]
  df_i  <- MPs_df[[mod_i]]
  
  for (j in module_names[-i]) {
    df_j <- MPs_df[[j]]
    shared <- intersect(df_i$genes, df_j$genes)
    
    for (g in shared) {
      freq_i <- df_i[g, "Freq"]
      freq_j <- df_j[g, "Freq"]
      
      if (freq_i > freq_j) {
        ## keep in i, remove from j
        MPs[[j]]     <- MPs[[j]][ MPs[[j]] != g ]
        MPs_df[[j]]  <- MPs_df[[j]][ rownames(MPs_df[[j]]) != g, , drop = FALSE ]
      } else if (freq_i < freq_j) {
        ## keep in j, remove from i
        MPs[[mod_i]]    <- MPs[[mod_i]][ MPs[[mod_i]] != g ]
        MPs_df[[mod_i]] <- MPs_df[[mod_i]][ rownames(MPs_df[[mod_i]]) != g, , drop = FALSE ]
      } else {
        ## equal frequency → keep the one that appears first in order
        if (which(module_names == mod_i) < which(module_names == j)) {
          MPs[[j]]     <- MPs[[j]][ MPs[[j]] != g ]
          MPs_df[[j]]  <- MPs_df[[j]][ rownames(MPs_df[[j]]) != g, , drop = FALSE ]
        } else {
          MPs[[mod_i]]    <- MPs[[mod_i]][ MPs[[mod_i]] != g ]
          MPs_df[[mod_i]] <- MPs_df[[mod_i]][ rownames(MPs_df[[mod_i]]) != g, , drop = FALSE ]
        }
      }
    }
  }
}

## ------------------------- FILTER MODULES ≥ 3 GENES ----------------------- ##
keep <- lengths(MPs) >= 3
MPs    <- MPs[keep]
MPs_df <- MPs_df[keep]

## ------------------------- SAVE OUTPUT ------------------------------------ ##
saveRDS(MPs,    file.path(OUT_DIR, paste0(CELLTYPE, ".rds")))
saveRDS(MPs_df, file.path(OUT_DIR, paste0(CELLTYPE, "_df.rds")))

message("Orthogonal genes saved to: ", OUT_DIR)





