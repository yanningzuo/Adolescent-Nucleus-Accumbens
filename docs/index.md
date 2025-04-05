# HRM-THC snRNA-seq Analysis Pipeline Documentation

This document provides a detailed overview of the analysis workflow for the HRM-THC snRNA-seq project. The pipeline is modularized into several steps, each implemented as an R script in the `scripts/` directory. The scripts have been directly converted from your existing code files and include detailed header comments.

---

## 1. Raw Data Processing

*Note:*  
Raw data processing is performed externally using Cell Ranger. The processed output (e.g., `filtered_feature_bc_matrix.h5`) should be placed in the `data/processed/` folder for further analysis in R.

---

## 2. Quality Control and Doublet Removal

**Script:** `02_quality_control_doublet_removal.R`

**Purpose:**  
- Load the merged gene expression matrix (from Cell Ranger) into a Seurat object.
- Calculate QC metrics such as RNA counts, gene counts, and percentages for mitochondrial and ribosomal genes.
- Apply filtering criteria to remove low-quality nuclei.
- Detect and remove doublets using DoubletFinder.

**Key Steps:**
- Compute metrics using `PercentageFeatureSet()`.
- Filter cells based on thresholds (e.g., 500 < nFeature_RNA < 6000, nCount_RNA between 1000 and 16000, <2% mitochondrial genes).
- Perform doublet detection with a pN value of 0.25 and an assumed doublet rate of 5%.
- Save the quality-controlled Seurat object.

---

## 3. SCTransform Clustering and Cell Type Annotation

**Script:** `03_sct_clustering_annotation.R`

**Purpose:**  
- Normalize the quality-controlled data using SCTransform, regressing out unwanted variation (e.g., percent.mt, nCount_RNA).
- Perform dimensionality reduction (PCA and UMAP) and cluster the cells.
- Annotate clusters using canonical markers (e.g., Snap25, Syt1 for neurons; Gja1 for astrocytes).

**Key Steps:**
- Apply SCTransform normalization.
- Run PCA and UMAP to reduce dimensionality.
- Cluster cells using `FindClusters()` with an appropriate resolution (e.g., 0.7).
- Optionally, rename clusters based on marker expression.
- Save the annotated Seurat object.

---

## 4. Differential Gene Expression (DEG) Analysis

**Script:** `04_DEG_analysis.R`

**Purpose:**  
- Identify differentially expressed genes (DEGs) comparing THC-treated vs. Vehicle and Het vs. WT within each cell cluster.
- Filter DEGs based on significance thresholds (e.g., FDR < 0.05, fold change thresholds).

**Key Steps:**
- Use Seurat’s `FindMarkers()` function for DEG detection.
- Filter results to retain significant DEGs (based on z-scores and adjusted p-values).
- Save DEG results for downstream pathway enrichment and further analysis.

---

## 5. Cell-to-Cell Communication Analysis

**Script:** `12_cell_to_cell_communication.R`

**Purpose:**  
- Analyze cell–cell interactions using LIANA.
- Convert LIANA's consensus resource to mouse orthologs.
- Run LIANA, aggregate the results, and filter significant interactions (aggregate rank ≤ 0.01).

**Key Steps:**
- Set up LIANA with a custom external resource.
- Execute `liana_wrap()` and `liana_aggregate()` to obtain interaction data.
- Visualize the interactions (e.g., with a dot plot).
- Save the LIANA analysis results.

---

## Running the Pipeline

Before running the pipeline, ensure that:
- You change directory names
- The R packages and software dependencies are installed.

Execute the scripts in the order listed above to perform a complete analysis of your snRNA-seq data.
