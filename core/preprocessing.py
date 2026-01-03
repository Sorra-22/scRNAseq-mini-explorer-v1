import anndata as ad
import scanpy as sc


def run_preprocessing(
    adata: ad.AnnData,
    n_hvgs: int = 2000,
    n_pcs: int = 50
) -> ad.AnnData:
    """
    Run minimal preprocessing for scRNA-seq analysis.

    Steps are executed ONLY if missing:
    - Normalization + log1p
    - Highly variable gene selection
    - PCA
    - Neighbors
    - UMAP

    Modifies adata IN PLACE and returns it.
    """

    # -------------------------
    # Normalization + log1p
    # -------------------------
    if "log1p" not in adata.uns:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

    # -------------------------
    # Highly variable genes
    # -------------------------
    if "highly_variable" not in adata.var.columns:
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=n_hvgs,
            subset=False
        )

    # -------------------------
    # Scaling (only for PCA)
    # -------------------------
    if "X_pca" not in adata.obsm:
        sc.pp.scale(adata, max_value=10)
        sc.tl.pca(adata, n_comps=n_pcs)

    # -------------------------
    # Neighbors
    # -------------------------
    if "neighbors" not in adata.uns:
        sc.pp.neighbors(adata, n_pcs=n_pcs)

    # -------------------------
    # UMAP
    # -------------------------
    if "X_umap" not in adata.obsm:
        sc.tl.umap(adata)

    return adata
