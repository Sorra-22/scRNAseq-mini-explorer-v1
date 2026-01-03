import streamlit as st
import scanpy as sc

from core.io import load_anndata
from core.inspect import inspect_anndata
from core.qc import compute_qc_metrics, filter_cells
from core.preprocessing import run_preprocessing
from core.clustering import run_leiden, rank_marker_genes

from ui.sidebar import (
    dataset_controls,
    qc_controls,
    preprocessing_controls,
    clustering_controls,
    visualization_controls,
)


def main():
    # ------------------------------------------------------------------
    # App configuration
    # ------------------------------------------------------------------
    st.set_page_config(
        page_title="scRNA-seq Mini-Explorer",
        layout="wide"
    )

    st.title("scRNA-seq Mini-Explorer")
    st.markdown(
        """
        Lightweight, local, Scanpy-based explorer for small-to-medium scRNA-seq datasets.
        Designed for CPU-only execution and limited memory environments.
        """
    )

    # ------------------------------------------------------------------
    # Sidebar: dataset selection
    # ------------------------------------------------------------------
    dataset_opts = dataset_controls()

    adata = load_anndata(
        uploaded_file=dataset_opts["uploaded_file"],
        use_example=dataset_opts["use_example"]
    )

    if adata is None:
        st.info("Please upload an `.h5ad` file or select the example dataset.")
        return

    # ------------------------------------------------------------------
    # Inspect dataset state
    # ------------------------------------------------------------------
    state = inspect_anndata(adata)

    # ------------------------------------------------------------------
    # Dataset summary
    # ------------------------------------------------------------------
    st.subheader("Dataset Summary")

    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Cells", adata.n_obs)
    with col2:
        st.metric("Genes", adata.n_vars)
    with col3:
        st.metric(
            "Embeddings",
            int(state["has_pca"]) + int(state["has_umap"])
        )

    with st.expander("AnnData details"):
        st.markdown("**.obs columns**")
        st.write(list(adata.obs.columns))

        st.markdown("**.var columns**")
        st.write(list(adata.var.columns))

        st.markdown("**Detected preprocessing state**")
        st.json(state)

    # ------------------------------------------------------------------
    # QC metrics + filtering
    # ------------------------------------------------------------------
    st.subheader("Quality Control")

    adata = compute_qc_metrics(adata)

    qc_opts = qc_controls(adata.n_obs)

    adata_filt, qc_stats = filter_cells(
        adata,
        min_genes=qc_opts["min_genes"],
        max_genes=qc_opts["max_genes"],
        max_pct_mt=qc_opts["max_pct_mt"]
    )

    st.write("**QC filtering summary**")
    st.json(qc_stats)

    if adata_filt.n_obs == 0:
        st.error("All cells were filtered out. Adjust QC thresholds.")
        return

    # ------------------------------------------------------------------
    # Preprocessing (reuse-first)
    # ------------------------------------------------------------------
    st.subheader("Preprocessing")

    prep_opts = preprocessing_controls(state)

    adata_filt = run_preprocessing(
        adata_filt,
        n_hvgs=prep_opts["n_hvgs"],
        n_pcs=prep_opts["n_pcs"]
    )

    st.success("Preprocessing complete (existing results reused where available).")

    # ------------------------------------------------------------------
    # Clustering
    # ------------------------------------------------------------------
    st.subheader("Clustering")

    cluster_opts = clustering_controls(state)

    adata_filt = run_leiden(
        adata_filt,
        resolution=cluster_opts["resolution"],
        force=cluster_opts["force_recluster"]
    )

    st.write(
        f"Leiden clusters detected: "
        f"{adata_filt.obs['leiden'].nunique()}"
    )

    if cluster_opts["compute_markers"]:
        with st.spinner("Computing marker genes..."):
            adata_filt = rank_marker_genes(adata_filt)
        st.success("Marker genes computed.")

    # ------------------------------------------------------------------
    # Visualization
    # ------------------------------------------------------------------
    st.subheader("Visualization")

    viz_opts = visualization_controls(
        obs_columns=list(adata_filt.obs.columns),
        gene_names=list(adata_filt.var_names)
    )

    plot_color = viz_opts["color_by"]
    if viz_opts["gene"] is not None:
        plot_color = viz_opts["gene"]

    fig = sc.pl.umap(
        adata_filt,
        color=plot_color,
        show=False,
        return_fig=True
    )

    st.pyplot(fig)

    # ------------------------------------------------------------------
    # Marker gene tables
    # ------------------------------------------------------------------
    if "rank_genes_groups" in adata_filt.uns:
        st.subheader("Marker Genes")

        clusters = adata_filt.obs["leiden"].unique().tolist()
        selected_cluster = st.selectbox(
            "Select cluster",
            options=clusters
        )

        markers_df = sc.get.rank_genes_groups_df(
            adata_filt,
            group=selected_cluster
        )

        st.dataframe(
            markers_df.head(20),
            use_container_width=True
        )

        st.download_button(
            label="Download marker genes (CSV)",
            data=markers_df.to_csv(index=False),
            file_name=f"markers_cluster_{selected_cluster}.csv",
            mime="text/csv"
        )

    # ------------------------------------------------------------------
    # Download filtered AnnData
    # ------------------------------------------------------------------
    st.subheader("Download")

    st.download_button(
        label="Download filtered AnnData (.h5ad)",
        data=adata_filt.write_h5ad,
        file_name="filtered_anndata.h5ad"
    )


if __name__ == "__main__":
    main()
