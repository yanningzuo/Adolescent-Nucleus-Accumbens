# snRNA-seq Analysis Workflow and Pipeline Documentation

This document provides a detailed overview of the analysis pipeline used for the HRM-THC snRNA-seq project. Below is a breakdown of each major processing and analysis step.

---

## 1. Raw Data Processing

Raw sequencing data is processed using Cell Ranger (v6.0.2) and aligned to the mm10-3.0.0 reference genome. Gene expression matrices from 13 samples are integrated using the `aggr` function (Cell Ranger v6.0.1), resulting in a combined dataset with a total of 82,251 nuclei.

---

## 2. Quality Control and Doublet Removal

The merged gene expression matrix is imported into Seurat (v4.9.9.9058). Quality control criteria include:
- **RNA count:** Between 1,000 and 16,000.
- **Gene count:** Between 500 and 6,000.
- **Mitochondrial gene percentage:** < 2%.
- **RPS and RPL gene percentages:** < 1% each.

Doublet detection is performed with DoubletFinder (v2.0.3) using custom parameters:
- Dimensionality reduction is carried out using the first 20 PCs.
- A pN value of 0.25 and an assumed doublet rate of 5% are applied.
- Visualization outputs (UMAP, feature plots) are generated for each sample.

*(See `scripts/02_quality_control_doublet_removal.R` for full code.)*

---

## 3. Single-Nucleus Transcriptome Processing, Clustering, and Annotation

Normalization is performed using `SCTransform()` with regression on mitochondrial percentage, sequencing depth, and ribosomal gene percentages. Clustering steps include:
- **PCA:** Run with the first 50 principal components.
- **Clustering:** Neighbors are identified using `FindNeighbors()`, and clusters are formed with `FindClusters()` at a resolution of 0.7.
- **Visualization:** UMAP is used for dimensionality reduction (`RunUMAP()`).
- **Cell Type Annotation:** Canonical markers (e.g., Snap25, Syt1 for neurons; Gja1 for astrocytes) are used to annotate clusters. Clusters with high mitochondrial content or mixed markers are removed.

*(Refer to `scripts/03_sct_clustering_annotation.R` for implementation details.)*

---

## 4. Differential Gene Expression (DEG) Analysis

Differential expression analysis is conducted using Seurat’s `FindMarkers()` function:
- **Comparisons:** Treatment (THC vs. Vehicle) and Genotype (Reln -/+ vs. WT) are compared within each cell cluster.
- **Thresholds:** DEGs are defined as those with an absolute average log2 fold change > 0.25 and an FDR-adjusted p-value < 0.05.
- **Output:** Lists of significant DEGs for each cluster are generated.

*(See `scripts/04_DEG_analysis.R` for code.)*

---

## 5. Pathway Annotation

Pathway enrichment is performed using the enrichR R package with databases:
- KEGG_2021_Human
- GO_Biological_Process_2021
- BioCarta_2016
- Reactome_2022

Enrichment is based on the deviation from the expected rank of terms, and p-values are adjusted using the Benjamini-Hochberg method.

*(Detailed code in `scripts/05_pathway_annotation.R`.)*

---

## 6. Weighted Gene Co-expression Network Analysis (WGCNA)

WGCNA is used to generate signed co-expression networks:
- **Normalization:** Raw counts are normalized and log-transformed.
- **Network Construction:** An adjacency matrix is generated based on biweight midcorrelation, followed by the calculation of a topological overlap measure (TOM).
- **Module Detection:** Gene dendrograms are constructed using hierarchical clustering, and modules are detected with dynamic tree cutting.
- **Module Merging:** Modules with highly correlated eigengenes (ME correlation > 0.85) are merged.

*(For full details, refer to `scripts/06_WGCNA_analysis.R`.)*

---

## 7. Trait Module Correlation

Pearson’s correlation is calculated between module eigengenes and traits (THC treatment and Reln genotype). Modules with significant correlations (p < 0.05) are highlighted.

*(Implemented in `scripts/07_trait_module_correlation.R`.)*

---

## 8. Associating WGCNA Modules with Human GWAS Traits

GWAS summary statistics for various traits (e.g., Cannabis Use Disorder, Opioid Use Disorder, Schizophrenia) are integrated:
- **Mapping:** SNPs are mapped to genes based on proximity and NAc tissue-specific eQTLs/sQTLs.
- **Enrichment:** Modules are tested for enrichment of GWAS-implicated genes using Marker Set Enrichment Analysis (MSEA) in Mergeomics.

*(See `scripts/08_GWAS_module_association.R`.)*

---

## 9. Weighted Key Driver Analysis (wKDA)

wKDA is performed using the Mergeomics package to identify key drivers in the gene regulatory networks:
- **Network Construction:** A brain-specific Bayesian network is used.
- **Key Driver Identification:** Genes with significant downstream connectivity are identified as key drivers.
- **Visualization:** Results are visualized with Cytoscape.

*(Refer to `scripts/09_wKDA.R`.)*

---

## 10. Trajectory Analysis

Pseudotime trajectories are inferred using Monocle3:
- **Dimensionality Reduction:** UMAP is applied using the first 30 PCs.
- **Graph Learning:** The `learn_graph()` function constructs the trajectory graph.
- **Pseudotime Assignment:** Cells are ordered along the trajectory with `order_cells()`.
- **Visualization:** Gene expression trends along pseudotime are visualized using GAM smoothers.

*(See `scripts/10_trajectory_analysis.R` for full details.)*

---

## 11. Conservation Analysis

To compare cell type conservation across species, t-statistics are calculated for each gene using SCTransformed data. Pearson’s correlations of these t-statistics are computed between datasets from mouse and human NAc, as well as juvenile mouse dorsal striatum.

*(Detailed implementation in `scripts/11_conservation_analysis.R`.)*

---

## 12. Cell-to-Cell Communication Analysis

Cell-cell interactions are identified using LIANA:
- **Resource Conversion:** LIANA’s Consensus resource is converted to murine symbols.
- **Analysis:** LIANA is run with a custom orthologous resource, and significant interactions (aggregate rank ≤ 0.01) are visualized.

*(See `scripts/12_cell_to_cell_communication.R`.)*

---

## 13. Gene Regulatory Network (GRN) Analysis

Transcription factor inference is performed using SCENIC:
- **GRN Construction:** Genie3 is used to infer initial regulatory relationships.
- **Regulon Definition:** Co-expression modules are refined into regulons.
- **Scoring:** Regulon activity is calculated and visualized using AUC scores.
- **Specificity:** Regulon specificity scores (RSS) are computed to identify cell-type-specific regulators.

*(Implemented in `scripts/13_GRN_analysis.R`.)*

---

## Running the Analysis

Before running the pipeline, ensure that:
- All raw data files are in the appropriate `data/raw/` directory.
- Required R packages and external software (Cell Ranger, Cytoscape, etc.) are installed.
- File paths in the scripts are correctly adjusted to match your local environment.

Execute the scripts sequentially as numbered in the **scripts/** directory to reproduce the full analysis pipeline.

---

*For any questions or further instructions, please consult the README.md or contact the project maintainer.*
