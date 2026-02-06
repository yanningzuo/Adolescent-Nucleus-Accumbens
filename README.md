# Adolescent-Nucleus-Accumbens

This repository provides a comprehensive pipeline for processing and analyzing single-nucleus RNA-sequencing (snRNA-seq) data from heterozygous reeler mice treated with vehicle and THC during adolescence. The analysis includes raw data processing, quality control, normalization, differential expression, gene co-expression network modeling, and cross-age conservation analysis.

## Requirements

* **Software:**
* [Cell Ranger](https://support.10xgenomics.com/) (v6.0.1 used for aggregation)
* R (version 4.x or above)


* **Key R Packages:**
* **Seurat (v4.x):** Core analysis framework
* **DoubletFinder:** Doublet detection
* **WGCNA:** Co-expression network construction
* **LIANA:** Ligand-receptor communication
* **SCENIC:** Gene regulatory network inference
* **Monocle3:** Trajectory inference
* **Mergeomics:** MSEA and Key Driver Analysis



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



## Getting Started

1. **Clone the Repository:**
```sh
git clone https://github.com/NumaanFormoli/Adolescent-Nucleus-Accumbens.git
cd Adolescent-Nucleus-Accumbens

```


2. **Prepare the Environment:**
* Install required R packages and ensure Cell Ranger is in your PATH.


3. **Run the Pipeline:**
* Execute each R script in sequential order (01 to 11).
* **Note:** Update `setwd()` and file paths in each script to match your local server environment.



## License

This project is licensed under the [MIT License](https://www.google.com/search?q=LICENSE).
