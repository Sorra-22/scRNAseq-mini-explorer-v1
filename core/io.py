from pathlib import Path
from typing import Optional

import streamlit as st
import scanpy as sc
import anndata as ad


# Path to example dataset (optional)
EXAMPLE_DATA_PATH = Path("data/pbmc3k_example.h5ad")


@st.cache_data(show_spinner=False)
def _load_h5ad_from_path(path: Path) -> ad.AnnData:
    """
    Load an AnnData object from a local path.

    Cached to avoid reloading on UI interactions.
    """
    return sc.read_h5ad(path)


@st.cache_data(show_spinner=False)
def _load_h5ad_from_upload(uploaded_file) -> ad.AnnData:
    """
    Load an AnnData object from a Streamlit uploaded file.
    """
    return sc.read_h5ad(uploaded_file)


def load_anndata(
    uploaded_file=None,
    use_example: bool = False
) -> Optional[ad.AnnData]:
    """
    Load AnnData based on user selection.

    Returns:
        AnnData or None if no dataset is selected.
    """
    if use_example:
        if not EXAMPLE_DATA_PATH.exists():
            st.error(
                "Example dataset not found. "
                "Please upload your own .h5ad file."
            )
            return None
        return _load_h5ad_from_path(EXAMPLE_DATA_PATH)

    if uploaded_file is not None:
        return _load_h5ad_from_upload(uploaded_file)

    return None


def save_anndata(adata: ad.AnnData, filename: str) -> Path:
    """
    Save AnnData to disk and return the path.

    Intended for filtered datasets.
    """
    output_path = Path(filename)
    adata.write(output_path)
    return output_path
