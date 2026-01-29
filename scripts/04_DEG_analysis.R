# Differentially expressed gene identification
library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)
library(patchwork)
set.seed(12)
# parallelization
library(future)
plan("multicore", workers = 8)
plan()
options(future.globals.maxSize = 100 * 1024 ^ 3) # for 100 Gb RAM
load("20230107_SCT_noGlut-61911.RData")
# remove ribosomal and mitochondrial genes
DefaultAssay(data) <- "SCT"
counts <- GetAssayData(object = data, slot = "counts", assay = "RNA")
counts <- counts[-which(startsWith(rownames(counts), "Rps")), ]
counts <- counts[-which(startsWith(rownames(counts), "Rpl")), ]
counts <- counts[-which(startsWith(rownames(counts), "mt-")), ]
data <- subset(data, features = rownames(counts))
Idents(object = data) <- data@meta.data[["detailed.celltype"]]
## DEGs
# logfc.threshold = 0
reln.markers <- list()
for (clu in levels(Idents(data))){
    reln.markers[[clu]] <- FindMarkers(object = data, logfc.threshold = 0, ident.1 = "Het", group.by = 'Genotype', subset.ident = clu)
}
names(reln.markers) <- levels(Idents(data))
thc.markers <- list()
for (clu in levels(Idents(data))){
    thc.markers[[clu]] <- FindMarkers(object = data, logfc.threshold = 0, ident.1 = "THC", group.by = 'Treatment', subset.ident = clu)
}
names(thc.markers) <- levels(Idents(data))
## significnace FDR < 5%, |z| > 1
# Bar plot Reln DEGs
sig.reln.markers <- list()
for (clu in levels(Idents(data))){
    temp <- reln.markers[[clu]]
    temp$z <- (temp$avg_log2FC - mean(temp$avg_log2FC))/sd(temp$avg_log2FC)
    temp <- temp[which(temp$p_val_adj < 0.05 & abs(temp$z) > 1), ]
    sig.reln.markers[[clu]] <- temp
}
names(sig.reln.markers) <- levels(Idents(data))
sig.thc.markers <- list()
for (clu in levels(Idents(data))){
    temp <- thc.markers[[clu]]
    temp$z <- (temp$avg_log2FC - mean(temp$avg_log2FC))/sd(temp$avg_log2FC)
    temp <- temp[which(temp$p_val_adj < 0.05 & abs(temp$z) > 1 ), ]
    sig.thc.markers[[clu]] <- temp
}
names(sig.thc.markers) <- levels(Idents(data))
