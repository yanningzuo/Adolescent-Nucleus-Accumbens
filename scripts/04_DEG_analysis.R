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

# Pathway analysis
library(enrichR)
dbs <- c("KEGG_2021_Human", "GO_Biological_Process_2021", "BioCarta_2016", "Reactome_2022")
reln.deg.go.up <- list()
for (clu in levels(Idents(data))){
    temp <- sig.reln.markers[[clu]]
    temp <- temp[which(temp$avg_log2FC > 0), ]
    genes <- rownames(temp)
    go.genes <- enrichr(genes, dbs)
    if (is.null(go.genes) == TRUE){
      reln.deg.go.up[[clu]] <- 0
    }else{
      list.go.genes = do.call("rbind",go.genes)
      list.go.genes$padj <- p.adjust(list.go.genes$P.value, method = "BH")
      reln.deg.go.up[[clu]] <- list.go.genes[order(list.go.genes$padj), ]
    }
}
reln.deg.go.down <- list()
for (clu in levels(Idents(data))){
    temp <- sig.reln.markers[[clu]]
    temp <- temp[which(temp$avg_log2FC < 0), ]
    genes <- rownames(temp)
    go.genes <- enrichr(genes, dbs)
    if (is.null(go.genes) == TRUE){
      reln.deg.go.down[[clu]] <- 0
    }else{
      list.go.genes = do.call("rbind", go.genes)
      list.go.genes$padj <- p.adjust(list.go.genes$P.value, method = "BH")
      reln.deg.go.down[[clu]] <- list.go.genes[order(list.go.genes$padj), ]
    }
}
thc.deg.go.up <- list()
for (clu in levels(Idents(data))){
    temp <- sig.thc.markers[[clu]]
    temp <- temp[which(temp$avg_log2FC > 0), ]
    genes <- rownames(temp)
    go.genes <- enrichr(genes, dbs)
    if (is.null(go.genes) == TRUE){
      thc.deg.go.up[[clu]] <- 0
    }else{
      list.go.genes = do.call("rbind", go.genes)
      list.go.genes$padj <- p.adjust(list.go.genes$P.value, method = "BH")
      thc.deg.go.up[[clu]] <- list.go.genes[order(list.go.genes$padj), ]
    }
}
thc.deg.go.down <- list()
for (clu in levels(Idents(data))){
    temp <- sig.thc.markers[[clu]]
    temp <- temp[which(temp$avg_log2FC < 0), ]
    genes <- rownames(temp)
    go.genes <- enrichr(genes, dbs)
    if (is.null(go.genes) == TRUE){
      thc.deg.go.down[[clu]] <- 0
    }else{
      list.go.genes = do.call("rbind", go.genes)
      list.go.genes$padj <- p.adjust(list.go.genes$P.value, method = "BH")
      thc.deg.go.down[[clu]] <- list.go.genes[order(list.go.genes$padj), ]
    }
}
save.image("DEG.RData")
