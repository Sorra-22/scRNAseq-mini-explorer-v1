import streamlit as st
from typing import Dict, Any, List


def dataset_controls() -> Dict[str, Any]:
    """
    Dataset selection controls.
    """
    st.sidebar.header("Dataset")

    source = st.sidebar.radio(
        "Data source",
        ["Upload .h5ad", "Example: PBMC 3k"]
    )

    uploaded_file = None
    use_example = False

    if source == "Upload .h5ad":
        uploaded_file = st.sidebar.file_uploader(
            "Upload AnnData file",
            type=["h5ad"]
        )
    else:
        use_example = True

    return {
        "uploaded_file": uploaded_file,
        "use_example": use_example
    }


def qc_controls(adata_n_obs: int) -> Dict[str, Any]:
    """
    QC filtering controls.
    """
    st.sidebar.header("QC Filtering")

    min_genes = st.sidebar.slider(
        "Min genes per cell",
        min_value=0,
        max_value=5000,
        value=200,
        step=50
    )

    max_genes = st.sidebar.slider(
        "Max genes per cell",
        min_value=500,
        max_value=10000,
        value=min(2500, 10000),
        step=100
    )

    max_pct_mt = st.sidebar.slider(
        "Max % mitochondrial",
        min_value=0.0,
        max_value=50.0,
        value=10.0,
        step=1.0
    )

    return {
        "min_genes": min_genes,
        "max_genes": max_genes,
        "max_pct_mt": max_pct_mt
    }


def preprocessing_controls(state: Dict[str, bool]) -> Dict[str, Any]:
    """
    Controls related to preprocessing reuse / execution.
    """
    st.sidebar.header("Preprocessing")

    reuse_pca = state.get("has_pca", False)
    reuse_umap = state.get("has_umap", False)

    st.sidebar.checkbox(
        "Reuse existing PCA",
        value=reuse_pca,
        disabled=True
    )

    st.sidebar.checkbox(
        "Reuse existing UMAP",
        value=reuse_umap,
        disabled=True
    )

    n_hvgs = st.sidebar.number_input(
        "Number of HVGs",
        min_value=500,
        max_value=5000,
        value=2000,
        step=500
    )

    n_pcs = st.sidebar.number_input(
        "Number of PCs",
        min_value=10,
        max_value=100,
        value=50,
        step=10
    )

    return {
        "n_hvgs": n_hvgs,
        "n_pcs": n_pcs
    }


def clustering_controls(state: Dict[str, bool]) -> Dict[str, Any]:
    """
    Clustering and marker gene controls.
    """
    st.sidebar.header("Clustering")

    resolution = st.sidebar.slider(
        "Leiden resolution",
        min_value=0.1,
        max_value=2.0,
        value=0.5,
        step=0.1
    )

    force_recluster = st.sidebar.checkbox(
        "Recompute Leiden clustering",
        value=False,
        disabled=not state.get("has_leiden", False)
    )

    compute_markers = st.sidebar.button(
        "Compute marker genes"
    )

    return {
        "resolution": resolution,
        "force_recluster": force_recluster,
        "compute_markers": compute_markers
    }


def visualization_controls(
    obs_columns: List[str],
    gene_names: List[str]
) -> Dict[str, Any]:
    """
    Visualization controls.
    """
    st.sidebar.header("Visualization")

    color_by = st.sidebar.selectbox(
        "Color UMAP by",
        options=obs_columns
    )

    gene = st.sidebar.selectbox(
        "Gene expression",
        options=["None"] + gene_names
    )

    return {
        "color_by": color_by,
        "gene": gene if gene != "None" else None
    }
