from typing import Optional
import anndata as ad
import scanpy as sc


def run_leiden(
    adata: ad.AnnData,
    resolution: float = 0.5,
    force: bool = False
) -> ad.AnnData:
    """
    Run Leiden clustering if not already present.

    Parameters:
        resolution: Leiden resolution parameter
        force: recompute even if Leiden already exists

    Modifies adata IN PLACE and returns it.
    """

    if "leiden" in adata.obs and not force:
        return adata

    sc.tl.leiden(
        adata,
        resolution=resolution,
        key_added="leiden"
    )

    return adata


def rank_marker_genes(
    adata: ad.AnnData,
    groupby: str = "leiden",
    method: str = "wilcoxon",
    force: bool = False
) -> Optional[ad.AnnData]:
    """
    Rank marker genes for clusters.

    Parameters:
        groupby: column in adata.obs
        method: statistical test
        force: recompute markers if already present

    Returns:
        AnnData or None if ranking cannot be performed
    """

    if groupby not in adata.obs:
        return None

    if "rank_genes_groups" in adata.uns and not force:
        return adata

    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        method=method
    )

    return adata
