import scanpy as sc
import os

adata = sc.read_10x_mtx(
    "filtered_gene_bc_matrices/hg19",
    var_names="gene_symbols",
    cache=True
)

adata.var_names_make_unique()
adata.write("pbmc3k.h5ad")

