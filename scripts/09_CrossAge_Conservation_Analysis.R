### conservation analyses

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

set.seed(12)

# parallelization
library(future)
plan("multicore", workers = 8)

## Conservation analysis between Adult NAc P70 and adolescent NAc

# Adolescent NAc t-stat prep
load("~/Desktop/HRM-THC-snRNA-seq/Merged/20230123_SCT_noGlut-61911.RData")
setwd("~/Desktop/HRM-THC-snRNA-seq/Merged")

DefaultAssay(data) <- "SCT"
Idents(object = data) <- data@meta.data$detailed.celltype

metadata <- FetchData(data, vars = "ident")
exp = data@assays[["SCT"]]@data

t.stats = data.frame(Genes = rownames(exp))
for (ident in unique(data@active.ident)){
    cells = rownames(metadata)[which(metadata$ident == ident)]
    temp.exp = exp[, which(colnames(exp) %in% cells)]
    other.cells = rownames(metadata)[which(metadata$ident != ident)]
    other.exp = exp[, which(colnames(exp) %in% other.cells)]

    s.celltype = apply(temp.exp, 1, sd)
    s.other = apply(other.exp, 1, sd)
    s.pooled = sqrt((s.celltype^2 + s.other^2)/2)
    m.celltype = apply(temp.exp, 1, mean)
    m.other = apply(other.exp, 1, mean)
    t = sqrt(length(cells)) * (m.celltype - m.other)/s.pooled

    t = as.data.frame(t)
    colnames(t) = ident
    t.stats = cbind(t.stats, t)

}

t.stats = t.stats[, -1]
save(t.stats, file = "t-stats.Rdata")

# Adult NAc P70 t-stat prep

setwd("~/Desktop/HRM-THC-snRNA-seq/Adult-NAc")

data <- read.csv("GSE118020_expression_matrix.csv.gz")
rownames(data) <- data$X
data <- data[, -1]
metadata <- read.csv("GSE118020_meta_data.csv.gz")

data <- CreateSeuratObject(counts = data, project = "Adult-NAc")
data <- AddMetaData(object = data, metadata = data.frame(detailed.ident = metadata$L2_Clustering, row.names = rownames(data@meta.data)))

data <- SCTransform(data, vars.to.regress = c("nCount_RNA"), verbose = TRUE)

DefaultAssay(data) <- "SCT"
Idents(object = data) <- data@meta.data[["detailed.ident"]]

metadata <- FetchData(data, vars = "ident")
exp = data@assays[["SCT"]]@data

t.stats = data.frame(Genes = rownames(exp))
for (ident in unique(data@active.ident)){
    cells = rownames(metadata)[which(metadata$ident == ident)]
    temp.exp = exp[, which(colnames(exp) %in% cells)]
    other.cells = rownames(metadata)[which(metadata$ident != ident)]
    other.exp = exp[, which(colnames(exp) %in% other.cells)]

    s.celltype = apply(temp.exp, 1, sd)
    s.other = apply(other.exp, 1, sd)
    s.pooled = sqrt((s.celltype^2 + s.other^2)/2)
    m.celltype = apply(temp.exp, 1, mean)
    m.other = apply(other.exp, 1, mean)
    t = sqrt(length(cells)) * (m.celltype - m.other)/s.pooled

    t = as.data.frame(t)
    colnames(t) = ident
    t.stats = cbind(t.stats, t)

}

t.stats = t.stats[, -1]
save(t.stats, file = "adult_t-stats.Rdata")

# conservation analysis
load("~/Desktop/HRM-THC-snRNA-seq/Merged/t-stats.RData")

adolescence = t.stats

load("~/Desktop/HRM-THC-snRNA-seq/Merged/adult_t-stats.RData")

adult = t.stats

shared.genes <- intersect(rownames(adolescence), rownames(adult))

adolescence <- adolescence[which(rownames(adolescence) %in% shared.genes), ]
adult <- adult[which(rownames(adult) %in% shared.genes), ]

# test with markers
load("~/Desktop/HRM-THC-snRNA-seq/Merged/20230123_SCT_noGlut-61911.RData")

marker = list()
for (ident in colnames(nac)){
    temp = markers[[ident]]
    temp = temp[which(temp$p_val_adj < 0.05 & temp$avg_log2FC > 0),]
    temp = temp[order(-temp$avg_log2FC), ]
    genes = rownames(temp)[1:100]
    marker[[ident]] = genes
}

marker = marker[-1]
marker = unique(marker)

adolescence <- adolescence[which(rownames(adolescence) %in% marker), ]
adult <- adult[which(rownames(adult) %in% marker), ]

cor <- ""
p <- ""
for (i in 1:dim(adolescence)[2]){
    temp.cor <- data.frame(cor = rep(0, dim(adult)[2]))
    temp.p <- data.frame(cor = rep(0, dim(adult)[2]))
    for (j in 1:dim(adult)[2]){
        cor.result = cor.test(adolescence[, i], adult[, j], method = c("pearson"))
        temp.cor[j, 1] = cor.result[["estimate"]]
        temp.p[j, 1] = cor.result[["p.value"]]
    }
    cor <- cbind(cor, temp.cor)
    p <- cbind(p, temp.p)
}

cor = cor[, -1]
p = p[, -1]

colnames(cor) = colnames(nac)
rownames(cor) = colnames(pfc)

colnames(p) = colnames(nac)
rownames(p) = colnames(pfc)

cor <- as.matrix(cor)
p <- as.matrix(p)

test_labels <- round(cor, 2)
test_labels[which(p >= 0.05)] = " "

cor <- as.data.frame(cor)

library(pheatmap)
library(RColorBrewer)
library(stringr)

breaksList = seq(-0.5, 0.5, by = 0.01)
pheatmap(cor, scale = "none", trace = "none", breaks = breaksList, col = colorRampPalette(rev(brewer.pal(9, "RdBu")))(length(breaksList)), cluster_cols = TRUE, cluster_rows = TRUE, border_color = NA, display_numbers = test_labels, fontsize_number = 20)
