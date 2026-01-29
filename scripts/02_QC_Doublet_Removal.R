setwd("~/Desktop/HRM-THC-snRNA-seq/Merged")
library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)
library(patchwork)
library(DoubletFinder)
library(here)
library(EnhancedVolcano)
set.seed(12)

# parallelization
library(future)
plan("multicore", workers = 8)
plan()
options(future.globals.maxSize = 100 * 1024 ^ 3) # for 100 Gb RAM

# read in data and stash metadata
counts <- Read10X_h5("filtered_feature_bc_matrix.h5")
metadata <- read.csv("aggregation.csv", stringsAsFactors=T)
data <- CreateSeuratObject(counts = counts, project = "HRM-THC") ## QC with 82251 cells before filtering

gemgroup <- sapply(strsplit(rownames(data@meta.data), split="-"), "[[", 2)
current.gemgroups <- seq(1, length(levels(metadata$sample_id)))
orig.ident <- levels(metadata$sample_id)
sampleID <- plyr::mapvalues(gemgroup, from = current.gemgroups, to = orig.ident)
data <- AddMetaData(object=data, metadata=data.frame(orig.ident=sampleID, row.names=rownames(data@meta.data)))
df <- data.frame(orig.ident=data@meta.data[["orig.ident"]], row.names=rownames(data@meta.data))
colnames(df) <- "Genotype"
df$Genotype[which(df$Genotype %in% metadata$sample_id[which(metadata$Genotype == "wt")])] <- "WT"
df$Genotype[which(df$Genotype %in% metadata$sample_id[which(metadata$Genotype == "het")])] <- "Het"
data <- AddMetaData(object=data, metadata=df)
df <- data.frame(orig.ident=data@meta.data[["orig.ident"]], row.names=rownames(data@meta.data))
colnames(df) <- "Treatment"
df$Treatment[which(df$Treatment %in% metadata$sample_id[which(metadata$Treatment == "VEH")])] <- "Veh"
df$Treatment[which(df$Treatment %in% metadata$sample_id[which(metadata$Treatment == "THC")])] <- "THC"
data <- AddMetaData(object=data, metadata=df)


data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^mt-")
data <- PercentageFeatureSet(data, pattern = "^Rpl", col.name = "percent.rpl")
data <- PercentageFeatureSet(data, pattern = "^Rps", col.name = "percent.rps")
VlnPlot(object = data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rps", "percent.rpl"), ncol = 5, group.by = "orig.ident", pt.size=0)

# filter the aggregated dataset for low quality cells, 72546 cells after filtering
data <- subset(data, subset = nFeature_RNA > 500
                  & nFeature_RNA < 6000
                  & nCount_RNA > 1000
                  & nCount_RNA < 16000
                  & percent.mt < 2
                  & percent.rps < 1
                  & percent.rpl < 1)

# Doublet removal with the assumption that doublets represent 5% of cells.
## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
FindDoublets <- function(library_id, seurat_aggregate) {
  rnaAggr <- seurat_aggregate
  seurat_obj <- subset(rnaAggr, idents = library_id)
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  seurat_obj <- RunPCA(seurat_obj)
  # ElbowPlot(seurat_obj)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
  DimPlot(seurat_obj)
  ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
  sweep.res.list <- paramSweep_v3(seurat_obj, PCs = 1:20, sct = F)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK <- bcmvn %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
    filter(BCmetric == max(BCmetric)) %>%
    select(pK)
  pK <- as.numeric(as.character(pK[[1]]))
  seurat_doublets <- doubletFinder_v3(seurat_obj, PCs = 1:20, pN = 0.25, pK = pK,
                               nExp = round(0.05*length(seurat_obj@active.ident)),
                               reuse.pANN = FALSE, sct = F)
  # create doublet groupings and visualize results
  DF.class <- names(seurat_doublets@meta.data) %>% str_subset("DF.classifications")
  pANN <- names(seurat_doublets@meta.data) %>% str_subset("pANN")
  p1 <- ggplot(bcmvn, aes(x=pK, y=BCmetric)) +
    geom_bar(stat = "identity") +
    ggtitle(paste0("pKmax=",pK)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p2 <- DimPlot(seurat_doublets, group.by = DF.class)
  p3 <- FeaturePlot(seurat_doublets, features = pANN)
  outFile <- paste0("doublets.",library_id,".pdf")
  dir.create("plots", showWarnings = FALSE)
  pdf(here("plots",outFile))
    print(p1) # need to use print() when drawing pdf in a function call
    print(p2)
    print(p3)
  dev.off()
  # create a df of barcodes and doublet designations
  df_doublet_barcodes <- as.data.frame(cbind(rownames(seurat_doublets@meta.data), seurat_doublets@meta.data[[DF.class]]))
  return(df_doublet_barcodes)
}

# take an aggregated snRNA seurat object and a list library_id and find doublets. return a df of doublet barcodes
# send DimPlot and FeaturePlot of doublets for each library to here("plots")
Idents(data) <- "orig.ident"
list.doublet.bc <- lapply(orig.ident, function(x) {FindDoublets(x, seurat_aggregate = data)})
doublet_id <- list.doublet.bc %>%
  bind_rows() %>%
  dplyr::rename("doublet_id" = "V2") %>%
  tibble::column_to_rownames(var = "V1") # this is the barcode column
table(doublet_id) # quantify total doublet vs. singlet calls (expect ~5% doublets)
# add doublet calls to aggregated snRNA object as doublet_id in meta.data slot
data <- AddMetaData(data, doublet_id)

# filter out doublets prior to snRNA preprocessing
Idents(data) <- "doublet_id"
data <- subset(data, idents = "Singlet")
# 68917 cells after QC filtering and doublet removal

save.image("20221110_68917.RData")
