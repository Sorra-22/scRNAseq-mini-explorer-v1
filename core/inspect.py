from typing import Dict
import anndata as ad


def inspect_anndata(adata: ad.AnnData) -> Dict[str, bool]:
    """
    Inspect an AnnData object and report the presence of
    common scRNA-seq preprocessing artifacts.

    This function MUST:
    - Never modify adata
    - Never compute anything
    - Be safe for any valid AnnData object
    """

    state = {}

    # -------------------------
    # Core structure
    # -------------------------
    state["has_X"] = adata.X is not None
    state["is_backed"] = adata.isbacked

    # -------------------------
    # QC metrics
    # -------------------------
    state["has_n_genes"] = "n_genes_by_counts" in adata.obs
    state["has_total_counts"] = "total_counts" in adata.obs
    state["has_pct_mt"] = "pct_counts_mt" in adata.obs

    # -------------------------
    # Gene-level annotations
    # -------------------------
    state["has_hvg"] = "highly_variable" in adata.var

    # -------------------------
    # Dimensionality reduction
    # -------------------------
    state["has_pca"] = "X_pca" in adata.obsm
    state["has_umap"] = "X_umap" in adata.obsm
    state["has_tsne"] = "X_tsne" in adata.obsm  # optional, not used

    # -------------------------
    # Graph-based analysis
    # -------------------------
    state["has_neighbors"] = "neighbors" in adata.uns

    # -------------------------
    # Clustering
    # -------------------------
    state["has_leiden"] = "leiden" in adata.obs
    state["has_louvain"] = "louvain" in adata.obs

    # -------------------------
    # Marker genes
    # -------------------------
    state["has_rank_genes"] = "rank_genes_groups" in adata.uns

    # -------------------------
    # Cell-type annotations (heuristic)
    # -------------------------
    state["has_celltype_annotations"] = any(
        col.lower() in ["cell_type", "celltype", "celltypes", "annotation", "annotations"]
        for col in adata.obs.columns
    )

    return state
