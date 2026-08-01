# Adolescent-Nucleus-Accumbens

This repository provides a comprehensive pipeline for processing and analyzing single-nucleus RNA-sequencing (snRNA-seq) data from heterozygous reeler mice treated with vehicle and THC during adolescence. The analysis includes raw data processing, quality control, normalization, differential expression, gene co-expression network modeling, cross-age conservation analysis, cell-cell communication analysis, and DEG analysis using bootstrapping.

## Requirements

**Software:**
* [Cell Ranger](https://support.10xgenomics.com/) (v6.0.1 used for aggregation)
* R (version 4.x or above)
* Python (version 3.10.18)

**Key R Packages:**
* **Seurat (v4.x):** Core analysis framework
* **DoubletFinder:** Doublet detection
* **WGCNA:** Co-expression network construction
* **LIANA:** Ligand-receptor communication
* **SCENIC:** Gene regulatory network inference
* **Monocle3:** Trajectory inference
* **Mergeomics:** MSEA and Key Driver Analysis
* **DESeq2:** pseudobulk DESeq analysis

**Key Python Packages:**
* **Scanpy (v. 1.9.3):** Python Single Cell analysis framework
* **Liana+ (v. 1.5.1):** Cell to Cell comunication inference
* **cell2cell (v. 0.7.4):** Tensor Cell2Cell factorization analysis
* **torch (v. >=1.12.1):** Torch backend used for factorization analysis

## Repository Structure

* **docs/**: `index.md` workflow guide.
* **scripts/**:
* **01_Raw_Processing.R**: CellRanger aggregation.
* **02_QC_Doublet_Removal.R**: Filtering and DoubletFinder.
* **03_SCT_Clustering_Annotation.R**: SCTransform and cell typing.
* **04_DEG_Pathway_Analysis.R**: DEGs and EnrichR (Merged Step 05).
* **05_WGCNA_Module_Enrichment.R**: Network building and module GO.
* **06_WGCNA_Trait_Correlation.R**: Module-trait associations.
* **07_Mergeomics_MSEA_KDA.R**: GWAS integration and hub gene identification.
* **08_Monocle3_Trajectory_Analysis.R**: Pseudotime inference.
* **09_CrossAge_Conservation_Analysis.R**: Adolescent vs. Adult comparison.
* **10_LIANA_Communication_Analysis.R**: Ligand-receptor signaling.
* **11_SCENIC_Regulon_Analysis.R**: Transcription factor networks.
* **12_DEG_Bootstrapping/**: DEG analysis for treatment groups using bootstrapping. (Directory)
* **13_LIANA+_C2C_Analysis.ipynb:** Liana+ and cell2cell based cell-cell communication factorization analysis.
* **14_DESeq2_pseudobulk.R:** DESeq2_pseudobulk analysis



## Getting Started

1. **Clone the Repository:**
```sh
git clone https://github.com/ftlabucsd/Adolescent-Nucleus-Accumbens.git
cd Adolescent-Nucleus-Accumbens

```


2. **Prepare the Environment:**
* Install required R packages and ensure Cell Ranger is in your PATH.
* Install required Python packages.

3. **Run the Pipeline:**
* Execute each R script in sequential order (01 to 13).
* **Note:** Update `setwd()` and file paths in each script to match your local server environment.


## License

This project is licensed under the [MIT License](https://en.wikipedia.org/wiki/MIT_License).
