# 🧬 scRNA-seq Mini-Explorer
*A lightweight, local-first single-cell RNA-seq exploration tool built with Scanpy and Streamlit*

## Overview

**scRNA-seq Mini-Explorer** is a general-purpose, lightweight Streamlit application for exploring single-cell RNA-seq datasets stored as AnnData (`.h5ad`) objects.

Inspired by **Seurat Explorer** and **Scanpy UI**, this tool is intentionally designed to:

- Run entirely on a local machine
- Be CPU-only
- Stay comfortably within ≤ 5 GB RAM
- Support most small-to-medium scRNA-seq datasets (up to ~50–100k cells)
- Avoid unnecessary recomputation by reusing existing analysis results

This project is ideal for students, bioinformaticians, and research labs who want a fast, local scRNA-seq exploration interface.

---

## Key Features

### 📂 Data Loading
- Upload any `.h5ad` file
- Optional built-in PBMC 3k example dataset
- Automatic inspection of dataset structure

### 🔍 Dataset Inspection
- Number of cells and genes
- Available `.obs` metadata fields
- Detection of existing PCA, neighbors, UMAP, and clustering

### 🧹 Quality Control (QC)
- Interactive cell filtering using:
  - Percentage of mitochondrial counts
  - Number of genes per cell
- QC metrics computed only if missing
- Download filtered AnnData object

### ⚙️ Preprocessing
- Highly variable gene selection
- Scaling
- PCA computation (reused if already present)

### 🧠 Clustering & Embeddings
- Neighborhood graph construction
- Leiden clustering
- UMAP visualization
- Automatic reuse of existing results

### 🔬 Biological Exploration
- Cluster-wise marker gene identification
- Interactive gene expression visualization
- Support for pre-annotated cell types (if present in `.obs`)

### 💾 Downloads
- Filtered `.h5ad` file
- Marker gene tables (CSV)

---

## Core Design Principles

- Local-first execution (no cloud, no GPU)
- Minimal recomputation
- Modular and readable codebase
- GitHub- and portfolio-ready
- Compatible with low-resource machines
