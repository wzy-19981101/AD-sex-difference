#!/usr/bin/env python3
"""
Single-cell RNA-seq preprocessing pipeline
==========================================
This script loads per-cell-type Loom files, merges them, performs QC,
normalization, highly-variable gene selection, Harmony batch correction,
Leiden clustering, UMAP visualization and marker-gene dot-plots.

Author:  wzy
Date:    2025
License: MIT
"""

import os
from typing import Dict, List

import scanpy as sc
import harmonypy 

# ------------------------- CONFIGURATION ------------------------- #

# Groups and cell types to process
GROUPS: List[str] = ["normal_F"]
CELL_TYPES: List[str] = [
    "astrocytes",
    "excitatory neurons",
    "inhibitory neurons",
    "microglia",
    "oligodendrocytes",
    "OPCs",
]

# Input / output paths
BASE_INPUT = "/home/wzy/zyx_nmf/NC/scoreinput"
OUTPUT_DIR = "/home/wzy/zyx_nmf/ecosystem/leiden_umap"
INTEGRATED_DIR = "/home/wzy/zyx_nmf/ecosystem/integrated_data"
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(INTEGRATED_DIR, exist_ok=True)

# Marker genes for dot-plot
MARKER_GENES: Dict[str, List[str]] = {
    "astrocytes": ["AQP4", "SLC1A2"],
    "excitatory neurons": ["NRNG", "CAMK2A", "SLC17A7"],
    "inhibitory neurons": ["GAD1", "GAD2"],
    "microglia": ["CSF1R", "CD74", "C3", "HLA-DRA", "CX3CR1", "C1QB"],
    "oligodendrocytes": ["MBP", "MOBP", "PLP1", "MOG"],
    "OPCs": ["VCAN", "PCDH15", "MEGF11", "PDGFRA", "SOX10"],
}

# Sample-type mapping
SAMPLE_TYPE_MAP = {"1": "AD", "2": "normal"}

# ------------------------- FUNCTIONS ------------------------- #


def load_and_concatenate(group: str) -> sc.AnnData:
    """
    Load per-cell-type Loom files and the external reference,
    then concatenate into a single AnnData object.
    """
    # Per-cell-type files
    adatas = [
        sc.read_loom(os.path.join(BASE_INPUT, f"{ct}_{group}.loom"))
        for ct in CELL_TYPES
    ]
    adata = sc.concat(adatas, join="outer")

    # External reference
    ex_path = os.path.join(BASE_INPUT, f"{group}_ROSMAP_ex.loom")
    adata_ex = sc.read_loom(ex_path)
    adata = sc.concat([adata_ex, adata], join="outer")
    return adata


def qc_and_normalize(adata: sc.AnnData) -> None:
    """
    Basic QC, normalization, HVG selection and scaling (in-place).
    """
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(
        adata, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes=5000
    )
    adata.raw = adata.copy()
    sc.pp.scale(adata, max_value=10)


def harmony_correction_and_cluster(adata: sc.AnnData) -> None:
    """
    Run Harmony batch correction on 'subject', compute neighbours,
    UMAP and Leiden clustering (in-place).
    """
    sc.external.pp.harmony_integrate(adata, key="subject")
    sc.pp.neighbors(adata, use_rep="X_pca_harmony", key_added="harmony_neighbors")
    sc.tl.umap(adata, neighbors_key="harmony_neighbors")
    sc.tl.leiden(
        adata,
        resolution=1,
        key_added="clusters",
        neighbors_key="harmony_neighbors",
        flavor="igraph",
        n_iterations=2,
        directed=False,
    )


def plot_umap(adata: sc.AnnData, group: str) -> None:
    """Save UMAP coloured by Leiden clusters."""
    sc.settings.figdir = OUTPUT_DIR
    sc.pl.umap(
        adata,
        color="clusters",
        frameon=False,
        title=f"{group} Leiden Clustering (UMAP)",
        palette="tab20",
        save=f"_{group}_clusters_1.pdf",
    )


def plot_dotplot(adata: sc.AnnData, group: str) -> None:
    """Filter marker genes for presence and save dot-plot."""
    present_markers = {
        ct: [g for g in genes if g in adata.var_names]
        for ct, genes in MARKER_GENES.items()
        if any(g in adata.var_names for g in genes)
    }
    sc.settings.figdir = OUTPUT_DIR
    sc.pl.dotplot(
        adata,
        present_markers,
        groupby="clusters",
        dendrogram=True,
        save=f"_{group}_Dotplot_markers.pdf",
    )


# ------------------------- MAIN PIPELINE ------------------------- #


def main() -> None:
    for group in GROUPS:
        print(f"Processing group: {group}")

        adata = load_and_concatenate(group)
        print(adata)
        print("Data loaded")

        qc_and_normalize(adata)
        print("QC & normalization finished")

        harmony_correction_and_cluster(adata)
        print("Harmony integration & clustering finished")

        # Map sample_type to human-readable labels
        adata.obs["sample_type"] = adata.obs["sample_type"].replace(SAMPLE_TYPE_MAP)

        # Save integrated object
        integrated_path = os.path.join(INTEGRATED_DIR, f"{group}.h5ad")
        adata.write(integrated_path)

        # Visualisations
        plot_umap(adata, group)
        plot_dotplot(adata, group)

        # Memory cleanup
        del adata


if __name__ == "__main__":
    main()
