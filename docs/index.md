# HRM-THC snRNA-seq Analysis Pipeline Documentation

This document provides a detailed overview of the analysis workflow for the HRM-THC snRNA-seq project. The pipeline is modularized into **14 sequential steps**, each implemented as an R script, directory, or notebook in the `scripts/` directory. 

---

## 1. Raw Data Processing

**Script:** `01_Raw_Processing.R`

* **Purpose:** Aggregates multiple 10x Genomics libraries into a single experiment using `cellranger aggr` without normalization to preserve biological variance for downstream SCTransform. 



## 2. Quality Control and Doublet Removal

**Script:** `02_QC_Doublet_Removal.R`

* **Purpose:** Filters the aggregated dataset for low-quality nuclei and removes doublets. 


* **Key Steps:** Filters based on `nFeature_RNA` (500-6000), `nCount_RNA` (1000-16000), and mitochondrial percentage (<2%). Doublets are identified and removed using the **DoubletFinder** package with a pN of 0.25. 



## 3. SCTransform Clustering and Annotation

**Script:** `03_SCT_Clustering_Annotation.R`

* **Purpose:** Normalizes data and defines cell identities. 


* **Key Steps:** Uses `SCTransform` to regress out mitochondrial and ribosomal bias. Clusters cells using 50 PCs and a resolution of 0.7 before annotating "coarse" and "detailed" cell types based on canonical markers like *Snap25*, *Gad1*, and *Drd1*. 



## 4. DEG and Pathway Analysis

**Script:** `04_DEG_Pathway_Analysis.R`

* **Purpose:** Identifies genes altered by Genotype or Treatment and interprets their biological meaning. 


* **Key Steps:** Runs `FindMarkers` to identify Differentially Expressed Genes (DEGs) for "Het" and "THC" conditions. Significant DEGs (FDR < 0.05 and ) are processed through **EnrichR** using KEGG, GO, and Reactome databases. 



## 5. WGCNA Module Enrichment

**Script:** `05_WGCNA_Module_Enrichment.R`

* **Purpose:** Constructs co-expression networks to find groups of highly correlated genes (modules). 


* **Key Steps:** Calculates soft-thresholding power for scale-free topology, builds a signed network, and identifies modules using dynamic tree cutting. Each module is then annotated for functional enrichment via **EnrichR**. 



## 6. WGCNA Trait Correlation

**Script:** `06_WGCNA_Trait_Correlation.R`

* **Purpose:** Statistically links gene modules to experimental metadata (Genotype/Treatment). 


* **Key Steps:** Calculates the correlation between module eigengenes (MEs) and external traits to identify which biological pathways are most affected by the experiment. 



## 7. Mergeomics (MSEA and KDA)

**Script:** `07_Mergeomics_MSEA_KDA.R`

* **Purpose:** Integrates WGCNA modules with external GWAS data (e.g., OUD). 


* **Key Steps:** Performs Marker Set Enrichment Analysis (**MSEA**) to find trait-associated modules and Key Driver Analysis (**KDA**) using a brain-specific network to identify "hub" genes. 



## 8. Monocle3 Trajectory Analysis

**Script:** `08_Monocle3_Trajectory_Analysis.R`

* **Purpose:** Models the developmental progression of interneurons, neuroblasts, and NPCs. 


* **Key Steps:** Infers a "pseudotime" path using **Monocle3**, allowing for the visualization of gene expression changes across a cellular differentiation gradient. 



## 9. Conservation Analysis (Cross-Age and Cross-Species)

**Script:** `09_Conservation_Analysis.R`

* **Purpose:** **Part 1 (cross-age):** compares the transcriptomic landscape of Adolescent NAc against Adult NAc (P70, GSE118020). **Part 2 (cross-species):** compares Adolescent Mouse NAc against Adult Human NAc (Tran et al., n = 8 donors) — manuscript Fig. S10.


* **Key Steps:** Generates cell-type specific t-statistics for each dataset and uses **Pearson Correlation** over the union of the top-100 upregulated marker genes per cell type to determine if cellular identities are conserved across development (Part 1) and across species via mouse-human ortholog mapping (Part 2). Both parts share the same statistical helper functions, defined once at the top of the script.



## 10. LIANA Communication Analysis

**Script:** `10_LIANA_Communication_Analysis.R`

* **Purpose:** Maps ligand-receptor signaling interactions between cell types. 


* **Key Steps:** Uses the **LIANA** framework to aggregate multiple interaction resources, specifically focusing on signaling between interneuron subpopulations and neuroblasts. 



## 11. SCENIC Regulon Analysis

**Script:** `11_SCENIC_Regulon_Analysis.R`

* **Purpose:** Reconstructs transcription factor (TF) regulatory networks. 


* **Key Steps:** Uses the **SCENIC** pipeline to identify "Regulons" (TFs and their targets) and calculates a **Regulon Specificity Score (RSS)** to identify TFs driving specific GABAergic identities. 

## 12. Bootstrapping based DEG Analysis
**Directory:** `12_DEG_Bootstrapping`
* **Purpose:** Find differentially-expressed genes between treatment and genotype groups.


* **Key Steps:** Samples 61,911 cells with replacement for 1,000 iterations such that sizes of groups was kept constant. Ran Seurat's FindMarkers between the groups to find DEGs.

## 13. Liana+ and cell2cell based cell-cell communication factorization analysis
**Notebook:** `13_LIANA+_C2C_Analysis.ipynb`
* **Purpose:** Determine how THC treatment affects intercellular communication patterns.


* **Key Steps:** Uses Liana+ to construct a 4 dimensional tensor across contexts (Treatments), sender cell types, reciever cell types, and ligand-receptor pairs.

## 14. DESeq2 Pseudobulk Analysis

**Script:** `14_DESeq2_pseudobulk.R`

* **Purpose:** Cross-validates the single-cell-level DEG calls (Section 4) with a pseudobulk approach, which is more robust to single-cell pseudoreplication.


* **Key Steps:** Aggregates raw counts per cell type per sample into pseudobulk profiles, then runs **DESeq2** (`~ libsize + Genotype + Treatment + Genotype:Treatment`) at both broad and detailed cell-type resolution to identify Treatment, Genotype, and Genotype x Treatment Interaction DEGs. Produces DEG count bar plots, volcano plots, and a hypergeometric overlap test against the single-cell DEG calls from Section 4.
---

## Running the Pipeline

Execute the scripts in the numerical order listed above (01 to 14). Ensure all file paths in the `load()` and `read.csv()` functions are updated to reflect your local environment or server directory structure.
