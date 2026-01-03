from typing import Tuple
import anndata as ad
import scanpy as sc


def compute_qc_metrics(adata: ad.AnnData) -> ad.AnnData:
    """
    Compute standard QC metrics if they are not already present.

    This function MAY modify adata in place but MUST NOT
    perform filtering or normalization.
    """

    obs = adata.obs

    needs_qc = not all(
        key in obs.columns
        for key in ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
    )

    if not needs_qc:
        return adata

    # Identify mitochondrial genes conservatively
    if "mt" not in adata.var.columns:
        adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt"],
        inplace=True
    )

    return adata


def filter_cells(
    adata: ad.AnnData,
    min_genes: int,
    max_genes: int,
    max_pct_mt: float
) -> Tuple[ad.AnnData, dict]:
    """
    Filter cells based on QC thresholds.

    Returns:
        filtered_adata: a COPY of AnnData
        stats: dictionary with filtering summary
    """

    initial_n = adata.n_obs

    mask = (
        (adata.obs["n_genes_by_counts"] >= min_genes) &
        (adata.obs["n_genes_by_counts"] <= max_genes) &
        (adata.obs["pct_counts_mt"] <= max_pct_mt)
    )

    filtered_adata = adata[mask].copy()

    stats = {
        "initial_cells": initial_n,
        "remaining_cells": filtered_adata.n_obs,
        "filtered_cells": initial_n - filtered_adata.n_obs
    }

    return filtered_adata, stats
