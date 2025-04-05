# Adolescent-Nucleus-Accumbens

# HRM-THC snRNA-seq Analysis Pipeline

This repository provides a comprehensive pipeline for processing and analyzing single-nucleus RNA-sequencing (snRNA-seq) data from HRM-THC experiments. The analysis includes raw data processing, quality control and doublet removal, normalization and clustering using SCTransform, differential gene expression (DEG) analysis, pathway enrichment, weighted gene co-expression network analysis (WGCNA), trait-module correlation, GWAS module association, weighted key driver analysis (wKDA), trajectory inference, conservation analysis, cell-to-cell communication, and gene regulatory network (GRN) analysis.

## Requirements

- **Software:**
  - [Cell Ranger](https://support.10xgenomics.com/) (v6.0.2 for raw processing, v6.0.1 for aggregation)
  - R (version 4.x or above)
  
- **R Packages:**
  - Seurat (v4.9.9.9058)
  - DoubletFinder (v2.0.3)
  - enrichR
  - WGCNA
  - Mergeomics
  - Monocle3
  - LIANA (v0.1.14)
  - SCENIC (v1.3.1)
  - plus other packages: dplyr, ggplot2, stringr, patchwork, future, here, etc.

*Note:* Please refer to each script’s header for additional package requirements and installation instructions.

## Repository Structure

- **data/**
  - **raw/**: Contains raw sequencing data and initial Cell Ranger outputs.
  - **processed/**: Contains merged gene expression matrices and quality-controlled data.

- **docs/**
  - **index.md**: Detailed documentation of the analysis workflow, including method descriptions and step-by-step pipeline instructions.

- **scripts/**
  - **01_raw_processing.R**: Processes raw sequencing data with Cell Ranger.
  - **02_quality_control_doublet_removal.R**: Performs quality control and doublet removal using Seurat and DoubletFinder.
  - **03_sct_clustering_annotation.R**: Conducts SCTransform normalization, PCA, clustering, and cell type annotation.
  - **04_DEG_analysis.R**: Identifies differentially expressed genes (DEGs) for treatment and genotype comparisons.
  - **05_pathway_annotation.R**: Performs pathway enrichment analysis using enrichR.
  - **06_WGCNA_analysis.R**: Builds weighted gene co-expression networks for each cell subtype.
  - **07_trait_module_correlation.R**: Correlates module eigengenes with treatment and genotype traits.
  - **08_GWAS_module_association.R**: Maps co-expression modules to human GWAS traits.
  - **09_wKDA.R**: Identifies key driver genes via weighted key driver analysis.
  - **10_trajectory_analysis.R**: Infers pseudotime trajectories using Monocle3.
  - **11_conservation_analysis.R**: Compares cell type profiles across species and studies.
  - **12_cell_to_cell_communication.R**: Analyzes cell-cell interactions using LIANA.
  - **13_GRN_analysis.R**: Constructs gene regulatory networks using SCENIC.

- **results/**
  - Stores all output files, figures, and tables from the analyses.

## Getting Started

1. **Clone the Repository:**
    ```sh
    git clone https://github.com/yourusername/HRM-THC-snRNAseq.git
    cd HRM-THC-snRNAseq
    ```

2. **Prepare the Environment:**
   - Install required software and R packages.
   - Place raw data in the `data/raw/` directory.

3. **Run the Pipeline:**
   - Execute each R script in sequential order as numbered in the `scripts/` directory.
   - Adjust file paths and parameters as needed based on your computing environment.

4. **Review Outputs:**
   - Analysis results will be saved in the `results/` folder.

For further details on each analysis step, please refer to the [Documentation](./docs/index.md).

## License

This project is licensed under the [MIT License](LICENSE).
