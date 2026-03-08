#!/usr/bin/env python
# -*- coding: utf-8 -*-
###############################################################################
##  09_modules_leiden_clustering.py
##
##  Perform Leiden clustering on NMF module weights.
##  1. Build AnnData from TXT files
##  2. PCA → neighbor → UMAP
##  3. Leiden clustering + filter cluster
##  4. Export UMAP plot & cluster cell lists
##
##  Author: wzy
##  Date  : 2025
##  Depends: scanpy, matplotlib, seaborn, bbknn, pyreadr
###############################################################################

import os
import glob
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import pyreadr
from bbknn import bbknn

## ------------------------------ PARAMETERS -------------------------------- ##
IN_DIR      = "/home/wzy/zyx_nmf/NC/score/microglia"   # *.txt files
OUT_DIR     = "/home/wzy/zyx_nmf/NC/cluster_programs"
CLUSTER_CSV = "microglia"                              # base name for outputs
FILTER_CL   = "7"                                      # cluster to remove (set None to skip)
sc.settings.figdir = OUT_DIR
sc.settings.set_figure_params(dpi=100, color_map='OrRd')

## --------------------- 1. LOAD NMF WEIGHTS → ANNDATA ---------------------- ##
txt_files = glob.glob(os.path.join(IN_DIR, "*.txt"))
data = {}
for f in txt_files:
    name = os.path.basename(f).replace(".txt", "")
    with open(f) as fh:
        genes = fh.read().strip().splitlines()
    data[name] = genes

df = pd.DataFrame.from_dict(data, orient='index').fillna('')
adata = sc.AnnData(df.notna().astype(float))   # 1/0 presence matrix
adata.obs_names = df.index

## parse obs meta
adata.obs['celltype_sample_disease_msex'] = adata.obs_names
adata.obs[['celltype', 'sample', 'disease', 'msex']] = (
    adata.obs['celltype_sample_disease_msex']
           .str.split('.', expand=True)
)

adata.raw = adata
adata.X = np.asarray(adata.X, dtype=float)

## --------------------- 2. PREPROCESS & DIMENSION REDUCTION ----------------- ##
sc.pp.highly_variable_genes(adata, subset=True)
sc.pp.pca(adata)
bbknn(adata, batch_key='sample')   # harmony alternative if preferred
sc.tl.umap(adata)

## --------------------- 3. LEIDEN CLUSTERING ------------------------------- ##
sc.tl.leiden(adata, resolution=1)

## plot initial
fig, ax = plt.subplots(figsize=(6, 5))
sc.pl.umap(adata, color='leiden', frameon=False,
           title='Leiden Clustering (all)', show=False, ax=ax)
fig.savefig(os.path.join(OUT_DIR, f'umap_{CLUSTER_CSV}_all.pdf'),
            bbox_inches='tight')

## --------------------- 4. FILTER CLUSTER (OPTIONAL) ----------------------- ##
if FILTER_CL is not None:
    adata_filt = adata[adata.obs['leiden'] != FILTER_CL].copy()
    sc.pp.neighbors(adata_filt, n_neighbors=50, use_rep='X_pca')
    sc.tl.umap(adata_filt)

    fig, ax = plt.subplots(figsize=(6, 5))
    sc.pl.umap(adata_filt, color='leiden', frameon=False,
               title=f'Leiden Clustering (cl {FILTER_CL} removed)',
               show=False, ax=ax)
    fig.savefig(os.path.join(OUT_DIR, f'umap_{CLUSTER_CSV}_filtered.pdf'),
                bbox_inches='tight')
else:
    adata_filt = adata

## --------------------- 5. EXPORT CLUSTER CELL LISTS ----------------------- ##
cluster_cells = {}
for cl in adata_filt.obs['leiden'].cat.categories:
    cells = adata_filt.obs_names[adata_filt.obs['leiden'] == cl].tolist()
    cluster_cells[cl] = cells

cluster_df = pd.DataFrame.from_dict(cluster_cells, orient='index').T
rds_path = os.path.join(OUT_DIR, f'{CLUSTER_CSV}.rds')
pyreadr.write_rds(rds_path, cluster_df)

print(f"Cluster cell lists saved to: {rds_path}")
