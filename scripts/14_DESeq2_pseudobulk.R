# =============================================================================
# DESeq2 Pseudobulk Differential Expression Analysis
# HRM-THC snRNA-seq Dataset
# Design: ~ libsize + Genotype + Treatment + Genotype:Treatment
# Contrasts: Treatment (THC vs Veh), Genotype (Het vs WT), Interaction
# =============================================================================


# =============================================================================
# Part 1: Pseudobulk data preparation
# =============================================================================

library(Seurat)
library(SingleCellExperiment)
library(matrixStats)
library(tidyverse)
library(knitr)
library(Matrix.utils)
library(magrittr)
library(DESeq2)

load("./Merged/20250324.RData")

## Extract Counts and Metadata

counts <- LayerData(data, layer = "counts", assay = "RNA")
metadata <- data@meta.data

Idents(object = data) <- "detailed.celltype"
# Idents(object = data) <- "coarse.celltype"
metadata$cluster_id <- factor(data@active.ident)

sce <- SingleCellExperiment(assays = list(counts = counts),
                           colData = metadata)

# Remove lowly expressed genes which have less than 10 cells with any counts
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]

## Preparing data for count aggregation

# Named vector of cluster names
kids <- purrr::set_names(levels(sce$cluster_id))

# Named vector of sample names
sids <- purrr::set_names(levels(as.factor(sce@colData@listData[["orig.ident"]])))

## Generate sample level metadata

# Determine the number of cells per sample
sce$sampleInfo <- sce@colData@listData[["orig.ident"]]

# Turn class "table" into a named vector of cells per sample
n_cells <- table(sce$sampleInfo) %>%  as.vector()
names(n_cells) <- names(table(sce$sampleInfo))

# Match the named vector with metadata to combine it
m <- match(names(n_cells), sce$sampleInfo)

# Create the sample level metadata by selecting specific columns
ei <- data.frame(colData(sce)[m, ],
                  n_cells, row.names = NULL) %>%
                dplyr::select("sampleInfo", "Genotype", "Treatment", "n_cells")

## Count aggregation
# Aggregate the counts per sample_id and cluster_id

# Subset metadata to only include the cluster and sample IDs to aggregate across
groups <- colData(sce)[, c("detailed.celltype", "sampleInfo")]
# groups <- colData(sce)[, c("coarse.celltype", "sampleInfo")]
groups$sampleInfo <- factor(groups$sampleInfo)

# Aggregate across cluster-sample groups
# Each row corresponds to aggregate counts for a cluster-sample combo
pb <- aggregate.Matrix(t(counts(sce)),
                       groupings = groups, fun = "sum")

# Not every cluster is present in all samples; create a vector that represents how to split samples
splitf <- sapply(stringr::str_split(rownames(pb),
                                    pattern = "_JB",n = 2), `[`, 1)

# Split data and turn into a list
# Each component corresponds to a cluster; storing associated expression matrix (counts)
# Transform data i.e, so rows are genes and columns are samples
pb <- split.data.frame(pb,factor(splitf)) %>%
        lapply(function(u)
                set_colnames(t(u), gsub(".*_", "", rownames(u))))

metadata <- metadata[, c(1, 4, 5, 16)]

save.image("./Merged/DEGs/20260504_DESeq2/20260504_pseudobulk-detailed.RData")

save.image("./Merged/DEGs/20260504_DESeq2/Broad/20260504_pseudobulk-broad.RData")


# =============================================================================
# Part 2: DESeq2 — Broad cell types (with reading depth as covariate)
# =============================================================================

library(DESeq2)

load("./Merged/DEGs/20260504_DESeq2/Broad/20260504_pseudobulk-broad.RData")
path = "./Merged/DEGs/20260504_DESeq2/Covariate_depth/Broad/"

for (cluster in names(pb)){
    work.path = paste0(file.path(path, cluster), "/")
    dir.create(work.path)
    setwd(work.path)

    # Subset the counts
    counts <- pb[[cluster]]
    counts@Dimnames[[2]] <- paste0("JB", counts@Dimnames[[2]])

    # Subset the metadata
    cluster_metadata <- metadata[which(metadata$cluster_id == cluster), -4]
    freq = as.data.frame(table(cluster_metadata$orig.ident))
    freq$Var1 <- gsub("_", "", freq$Var1)

    cluster_metadata <- cluster_metadata[!duplicated(cluster_metadata), ]
    rownames(cluster_metadata) <- cluster_metadata$orig.ident
    rownames(cluster_metadata) <- gsub("_", "", rownames(cluster_metadata))
    cluster_metadata <- cluster_metadata[, -1]
    cluster_metadata[["Treatment"]] <- factor(cluster_metadata[["Treatment"]], levels = c("Veh", "THC"))
    cluster_metadata[["Genotype"]] <- factor(cluster_metadata[["Genotype"]],  levels = c("WT", "Het"))
    cluster_metadata$n_cells = freq$Freq[match(freq$Var1, rownames(cluster_metadata))]
    cluster_metadata$n_cells = as.numeric(cluster_metadata$n_cells)
    cluster_metadata$libsize <- as.numeric(colSums(counts))
    cluster_metadata$libsize = scale(cluster_metadata$libsize)
    colnames(cluster_metadata)[4] = "libsize"

    # Check that all of the row names of the metadata are the same and in the same order as the column names of the counts in order to use as input to DESeq2
    cluster_counts <- data.frame(as.matrix(counts[, which(colnames(counts) %in% rownames(cluster_metadata))]))
    all(rownames(cluster_metadata) == colnames(cluster_counts))
    featureData <- data.frame(gene = rownames(cluster_counts))
    cluster_counts <- cbind(featureData, cluster_counts)

    dds <- DESeqDataSetFromMatrix(countData = cluster_counts,
                              colData = cluster_metadata,
                              design= ~ libsize + Genotype + Treatment + Genotype:Treatment, tidy = TRUE)

    # filter out genes where there are less than 3 samples with normalized counts greater than or equal to 5.

    dds <- estimateSizeFactors(dds)
    idx <- rowSums(counts(dds, normalized=TRUE) >= 5 ) >= 3
    dds <- dds[idx,]
    dds <- DESeq(dds)

    # Treatment effect
    res.thc <- results(dds, name = "Treatment_THC_vs_Veh")
    thc <- res.thc@listData
    thc <- data.frame(matrix(unlist(thc), nrow = length(thc), byrow = T), stringsAsFactors = FALSE)
    thc <- t(thc)
    rownames(thc) <- res.thc@rownames
    colnames(thc)<- c("basemean", "log2FC", "lfcSE", "stat", "p", "padj")
    thc <- data.frame(thc, stringsAsFactors = FALSE)
    thc <- thc[order(thc$padj), ]
    write.csv(thc, file = paste0(cluster, "_treatment-DEGs.csv"))

    # Genotype effect
    res.geno <- results(dds, name = "Genotype_Het_vs_WT")
    geno <- res.geno@listData
    geno <- data.frame(matrix(unlist(geno), nrow = length(geno), byrow = T), stringsAsFactors = FALSE)
    geno <- t(geno)
    rownames(geno) <- res.geno@rownames
    colnames(geno)<- c("basemean", "log2FC", "lfcSE", "stat", "p", "padj")
    geno <- data.frame(geno, stringsAsFactors = FALSE)
    geno <- geno[order(geno$padj), ]
    write.csv(geno, file = paste0(cluster, "_genotype-DEGs.csv"))

    # Treatment x Genotype effect
    res.int <- results(dds, name = "GenotypeHet.TreatmentTHC")
    int <- res.int@listData
    int <- data.frame(matrix(unlist(int), nrow = length(int), byrow = T), stringsAsFactors = FALSE)
    int <- t(int)
    rownames(int) <- res.int@rownames
    colnames(int)<- c("basemean", "log2FC", "lfcSE", "stat", "p", "padj")
    int <- data.frame(int, stringsAsFactors = FALSE)
    int <- int[order(int$padj), ]
    write.csv(int, file = paste0(cluster, "_int-DEGs.csv"))
}

## Compile results into single excel files

library(WriteXLS)

input_paths = "./Merged/DEGs/20260504_DESeq2/Covariate_depth/Broad/"
cell.types = list.files(input_paths)

result = list()
for (clu in cell.types){

    setwd(paste0(input_paths, clu, "/"))
    files = list.files(path = ".", pattern = ".csv")

    for(f in files){
        filename = gsub("-DEGs.csv", "", f)
        result[[filename]] <- read.csv(f)
    }
}

WriteXLS(result, "./Merged/DEGs/20260504_DESeq2/Covariate_depth/DESeq2_broad-cell-type_DEGs.xlsx", names(result), row.names = TRUE)


# =============================================================================
# Part 3: DESeq2 — Detailed cell types (with reading depth as covariate)
# =============================================================================

load("./Merged/DEGs/20260504_DESeq2/Detailed/20260504_pseudobulk-detailed.RData")
path = "./Merged/DEGs/20260504_DESeq2/Covariate_depth/"

for (cluster in names(pb)){
    work.path = paste0(file.path(path, cluster), "/")
    dir.create(work.path)
    setwd(work.path)

    # Subset the counts
    counts <- pb[[cluster]]
    counts@Dimnames[[2]] <- paste0("JB", counts@Dimnames[[2]])

    # Subset the metadata
    cluster_metadata <- metadata[which(metadata$cluster_id == cluster), -4]
    freq = as.data.frame(table(cluster_metadata$orig.ident))
    freq$Var1 <- gsub("_", "", freq$Var1)

    cluster_metadata <- cluster_metadata[!duplicated(cluster_metadata), ]
    rownames(cluster_metadata) <- cluster_metadata$orig.ident
    rownames(cluster_metadata) <- gsub("_", "", rownames(cluster_metadata))
    cluster_metadata <- cluster_metadata[, -1]
    cluster_metadata[["Treatment"]] <- factor(cluster_metadata[["Treatment"]], levels = c("Veh", "THC"))
    cluster_metadata[["Genotype"]] <- factor(cluster_metadata[["Genotype"]],  levels = c("WT", "Het"))
    cluster_metadata$n_cells = freq$Freq[match(freq$Var1, rownames(cluster_metadata))]
    cluster_metadata$n_cells = as.numeric(cluster_metadata$n_cells)
    cluster_metadata$libsize <- as.numeric(colSums(counts))
    cluster_metadata$libsize = scale(cluster_metadata$libsize)
    colnames(cluster_metadata)[4] = "libsize"

    # Check that all of the row names of the metadata are the same and in the same order as the column names of the counts in order to use as input to DESeq2
    cluster_counts <- data.frame(as.matrix(counts[, which(colnames(counts) %in% rownames(cluster_metadata))]))
    all(rownames(cluster_metadata) == colnames(cluster_counts))
    featureData <- data.frame(gene = rownames(cluster_counts))
    cluster_counts <- cbind(featureData, cluster_counts)

    dds <- DESeqDataSetFromMatrix(countData = cluster_counts,
                              colData = cluster_metadata,
                              design= ~ libsize + Genotype + Treatment + Genotype:Treatment, tidy = TRUE)

    # filter out genes where there are less than 3 samples with normalized counts greater than or equal to 5.

    dds <- estimateSizeFactors(dds)
    idx <- rowSums(counts(dds, normalized=TRUE) >= 5 ) >= 3
    dds <- dds[idx,]
    dds <- DESeq(dds)

    # Treatment effect
    res.thc <- results(dds, name = "Treatment_THC_vs_Veh")
    thc <- res.thc@listData
    thc <- data.frame(matrix(unlist(thc), nrow = length(thc), byrow = T), stringsAsFactors = FALSE)
    thc <- t(thc)
    rownames(thc) <- res.thc@rownames
    colnames(thc)<- c("basemean", "log2FC", "lfcSE", "stat", "p", "padj")
    thc <- data.frame(thc, stringsAsFactors = FALSE)
    thc <- thc[order(thc$padj), ]
    write.csv(thc, file = paste0(cluster, "_treatment-DEGs.csv"))

    # Genotype effect
    res.geno <- results(dds, name = "Genotype_Het_vs_WT")
    geno <- res.geno@listData
    geno <- data.frame(matrix(unlist(geno), nrow = length(geno), byrow = T), stringsAsFactors = FALSE)
    geno <- t(geno)
    rownames(geno) <- res.geno@rownames
    colnames(geno)<- c("basemean", "log2FC", "lfcSE", "stat", "p", "padj")
    geno <- data.frame(geno, stringsAsFactors = FALSE)
    geno <- geno[order(geno$padj), ]
    write.csv(geno, file = paste0(cluster, "_genotype-DEGs.csv"))

    # Treatment x Genotype effect
    res.int <- results(dds, name = "GenotypeHet.TreatmentTHC")
    int <- res.int@listData
    int <- data.frame(matrix(unlist(int), nrow = length(int), byrow = T), stringsAsFactors = FALSE)
    int <- t(int)
    rownames(int) <- res.int@rownames
    colnames(int)<- c("basemean", "log2FC", "lfcSE", "stat", "p", "padj")
    int <- data.frame(int, stringsAsFactors = FALSE)
    int <- int[order(int$padj), ]
    write.csv(int, file = paste0(cluster, "_int-DEGs.csv"))
}

## Compile results into single excel files

input_paths = "./Merged/DEGs/20260504_DESeq2/Covariate_depth/"
cell.types = list.files(input_paths)

result = list()
for (clu in cell.types){

    setwd(paste0(input_paths, clu, "/"))
    files = list.files(path = ".", pattern = ".csv")

    for(f in files){
        filename = gsub("-DEGs.csv", "", f)
        result[[filename]] <- read.csv(f)
    }
}

WriteXLS(result, "./Merged/DEGs/20260504_DESeq2/Covariate_depth/DESeq2_detailed-cell-type_DEGs.xlsx", names(result), row.names = TRUE)


# =============================================================================
# Part 4: Visualization
# =============================================================================

library(ggplot2)
library(patchwork)
library(dplyr)
library(WriteXLS)
library(ggrepel)

### DESeq2 with covariate (reads) pseudobulk result plots

## DEG list

path = "./Merged/DEGs/20260504_DESeq2/Covariate_depth/Detailed/"
cell.types = list.files(path)

thc.deg <- list()
for (clu in cell.types){

    setwd(paste0(path, clu, "/"))
    file = list.files(path = ".", pattern = "_treatment-DEGs.csv")
    temp <- read.csv(file)
    temp <- temp[which(temp$padj < 0.05), ]
    thc.deg[[clu]] <- temp

}

geno.deg <- list()
for (clu in cell.types){

    setwd(paste0(path, clu, "/"))
    file = list.files(path = ".", pattern = "_genotype-DEGs.csv")
    temp <- read.csv(file)
    temp <- temp[which(temp$padj < 0.05), ]
    geno.deg[[clu]] <- temp

}

int.deg <- list()
for (clu in cell.types){

    setwd(paste0(path, clu, "/"))
    file = list.files(path = ".", pattern = "_int-DEGs.csv")
    temp <- read.csv(file)
    temp <- temp[which(temp$padj < 0.05), ]
    int.deg[[clu]] <- temp

}

## DEG count bar plot

deg.count = data.frame(
  cell.type = names(thc.deg),
  treatment.up = rep("", length(names(thc.deg))),
  treatment.down = rep("", length(names(thc.deg))),
  genotype.up = rep("", length(names(thc.deg))),
  genotype.down = rep("", length(names(thc.deg))),
  int.up = rep("", length(names(thc.deg))),
  int.down = rep("", length(names(thc.deg)))
)

for (clu in cell.types){

    i = which(deg.count$cell.type == clu)
    deg.count$treatment.up[i] = dim(thc.deg[[clu]][which(thc.deg[[clu]]$log2FC > 0), ])[1]
    deg.count$treatment.down[i] = dim(thc.deg[[clu]][which(thc.deg[[clu]]$log2FC < 0), ])[1]
    deg.count$genotype.up[i] = dim(geno.deg[[clu]][which(geno.deg[[clu]]$log2FC > 0), ])[1]
    deg.count$genotype.down[i] = dim(geno.deg[[clu]][which(geno.deg[[clu]]$log2FC < 0), ])[1]
    deg.count$int.up[i] = dim(int.deg[[clu]][which(int.deg[[clu]]$log2FC > 0), ])[1]
    deg.count$int.down[i] = dim(int.deg[[clu]][which(int.deg[[clu]]$log2FC < 0), ])[1]

}

# stacked based on DEG types

df.up <- data.frame(
  group = rep(deg.count$cell.type, 3),
  layer = rep(c("Treatment", "Genotype", "Interaction"), each = length(deg.count$cell.type)),
  value = c(deg.count$treatment.up, deg.count$genotype.up, deg.count$int.up)
)

df.up$layer <- factor(df.up$layer, levels = c("Treatment", "Genotype", "Interaction"))
df.up$value <- as.numeric(df.up$value)

p1 <- ggplot(df.up, aes(x = group, y = value, fill = layer)) +
        geom_bar(stat = "identity", width = 0.6) +
        labs(
          title  = "Up-regulated DEGs",
          x = NULL,
          y = "Counts",
          fill = "DEG"
         ) +
        theme_minimal(base_size = 13) +
        theme(text = element_text(size = 16),
          axis.text.x = element_text(size = 12, angle = 45, hjust=1),
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
        )

df.down <- data.frame(
  group = rep(deg.count$cell.type, 3),
  layer = rep(c("Treatment", "Genotype", "Interaction"), each = length(deg.count$cell.type)),
  value = c(deg.count$treatment.down, deg.count$genotype.down, deg.count$int.down)
)

df.down$layer <- factor(df.down$layer, levels = c("Treatment", "Genotype", "Interaction"))
df.down$value <- as.numeric(df.down$value)

p2 <- ggplot(df.down, aes(x = group, y = value, fill = layer)) +
        geom_bar(stat = "identity", width = 0.6) +
        labs(
          caption  = "Down-regulated DEGs",
          x = NULL,
          y = "Counts",
          fill = "DEG"
         ) +
        theme_minimal(base_size = 13) +
        theme(text = element_text(size = 16),
          axis.text.x = element_text(size = 12, angle = 45, hjust=1),
          plot.caption = element_text(hjust = 0.5, size = 14, face = "bold")
        ) +
        scale_y_reverse() +
        ylim(300, 0)

p1 / p2

# stacked based on up/down

df.thc <- data.frame(
  group = rep(deg.count$cell.type, 2),
  layer = rep(c("Up", "Down"), each = length(deg.count$cell.type)),
  value = c(deg.count$treatment.up, deg.count$treatment.down)
)

df.thc$layer <- factor(df.thc$layer, levels = c("Up", "Down"))
df.thc$value <- as.numeric(df.thc$value)

p1 <- ggplot(df.thc, aes(x = group, y = value, fill = layer)) +
        geom_bar(stat = "identity", width = 0.6) +
        labs(
          title  = "Treatment DEGs",
          x = NULL,
          y = "Counts",
          fill = "Direction"
         ) +
        theme_minimal(base_size = 13) +
        theme(text = element_text(size = 16),
          axis.text.x = element_text(size = 12, angle = 45, hjust=1),
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
        )

df.geno <- data.frame(
  group = rep(deg.count$cell.type, 2),
  layer = rep(c("Up", "Down"), each = length(deg.count$cell.type)),
  value = c(deg.count$genotype.up, deg.count$genotype.down)
)

df.geno$layer <- factor(df.geno$layer, levels = c("Up", "Down"))
df.geno$value <- as.numeric(df.geno$value)

p2 <- ggplot(df.geno, aes(x = group, y = value, fill = layer)) +
        geom_bar(stat = "identity", width = 0.6) +
        labs(
          title  = "Genotype DEGs",
          x = NULL,
          y = "Counts",
          fill = "Direction"
         ) +
        theme_minimal(base_size = 13) +
        theme(text = element_text(size = 16),
          axis.text.x = element_text(size = 12, angle = 45, hjust=1),
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
        )

df.int <- data.frame(
  group = rep(deg.count$cell.type, 2),
  layer = rep(c("Up", "Down"), each = length(deg.count$cell.type)),
  value = c(deg.count$int.up, deg.count$int.down)
)

df.int$layer <- factor(df.int$layer, levels = c("Up", "Down"))
df.int$value <- as.numeric(df.int$value)

p3 <- ggplot(df.int, aes(x = group, y = value, fill = layer)) +
        geom_bar(stat = "identity", width = 0.6) +
        labs(
          title  = "Interaction DEGs",
          x = NULL,
          y = "Counts",
          fill = "Direction"
         ) +
        theme_minimal(base_size = 13) +
        theme(text = element_text(size = 16),
          axis.text.x = element_text(size = 12, angle = 45, hjust=1),
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
        )
p1 / p2 / p3

## volcano plot

# all DEGs - detailed cell type

path = "./Merged/DEGs/20260504_DESeq2/Covariate_depth/Detailed/"
cell.types = list.files(path)

thc.deg.all <- list()
for (clu in cell.types){

    setwd(paste0(path, clu, "/"))
    file = list.files(path = ".", pattern = "_treatment-DEGs.csv")
    thc.deg.all[[clu]] <- read.csv(file)

}

geno.deg.all <- list()
for (clu in cell.types){

    setwd(paste0(path, clu, "/"))
    file = list.files(path = ".", pattern = "_genotype-DEGs.csv")
    geno.deg.all[[clu]] <- read.csv(file)

}

int.deg.all <- list()
for (clu in cell.types){

    setwd(paste0(path, clu, "/"))
    file = list.files(path = ".", pattern = "_int-DEGs.csv")
    int.deg.all[[clu]] <- read.csv(file)

}

deg.all <- list()
for (clu in cell.types){

    thc = thc.deg.all[[clu]]
    thc$Type = rep("Treatment", dim(thc)[1])
    thc$Type[which(thc$padj >= 0.05)] = "NS"
    thc <- thc %>% distinct(Type, X, .keep_all = TRUE)
    geno = geno.deg.all[[clu]]
    geno$Type = rep("Genotype", dim(geno)[1])
    geno$Type[which(geno$padj >= 0.05)] = "NS"
    geno <- geno %>% distinct(Type, X, .keep_all = TRUE)
    int= int.deg.all[[clu]]
    int$Type = rep("Interaction", dim(int)[1])
    int$Type[which(int$padj >= 0.05)] = "NS"
    int <- int %>% distinct(Type, X, .keep_all = TRUE)
    deg.all[[clu]] <- rbind(thc, geno, int)

}

color_map <- c(
  Treatment   = "#E64B35",
  Genotype    = "#4DBBD5",
  Interaction = "#00A087",
  NS          = "grey75"
)

pics <- list()
for (clu in cell.types){

    df = deg.all[[clu]]
    col_custom <- setNames(color_map[df$Type], df$X)
    df$Type <- factor(df$Type,
                          levels = c("Treatment", "Genotype", "Interaction", "NS"))
    label_df <- df |>
      filter(Type != "NS") |>
      group_by(Type) |>
      slice_min(order_by = padj, n = 10) |>
      ungroup()

     p <- ggplot(df, aes(x = log2FC, y = -log10(padj), color = Type)) +
            geom_point(data = filter(df, Type == "NS"),
                        size = 1, alpha = 0.4) +
             geom_point(data = filter(df, Type != "NS"),
                        size = 1.5, alpha = 0.8) +
             geom_hline(yintercept = -log10(0.05),
                        linetype = "dashed", color = "grey40", linewidth = 0.4) +
             geom_vline(xintercept = c(-0.25, 0.25),
                        linetype = "dashed", color = "grey40", linewidth = 0.4) +
             geom_text_repel(
               data          = label_df,
               aes(label     = X),
               size          = 2.8,
               max.overlaps  = 20,
               box.padding   = 0.3,
               point.padding = 0.2,
               segment.color = "grey50",
               segment.size  = 0.3,
               show.legend   = FALSE
             ) +
             scale_color_manual(
               values = color_map,
               labels = c(
                 treatment   = "Treatment DEGs",
                 genotype    = "Genotype DEGs",
                 interaction = "Interaction DEGs",
                 NS          = "Not significant"
              )
             ) +
             labs(
               title = clu,
               x     = expression(log[2]~"Fold Change"),
               y     = expression(-log[10]~"(p"[adj]*")"),
               color = NULL
            ) +
             theme_classic(base_size = 13) +
             theme(
               legend.position  = "top",
               legend.key.size  = unit(0.4, "cm"),
               plot.title       = element_text(hjust = 0.5)
             ) +
             guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

    pics[[clu]] <- p

}

p <- Reduce( `+`, pics ) +
        patchwork::plot_layout(ncol = 5)
ggplot2::ggsave("20260529_volcano-DESeq2.pdf", plot = p, width = 20, height = 25, units = "in")

# all DEGs - broad cell type

path = "./Merged/DEGs/20260504_DESeq2/Covariate_depth/Broad/"
cell.types = list.files(path)

thc.deg.all <- list()
for (clu in cell.types){

    setwd(paste0(path, clu, "/"))
    file = list.files(path = ".", pattern = "_treatment-DEGs.csv")
    thc.deg.all[[clu]] <- read.csv(file)

}

geno.deg.all <- list()
for (clu in cell.types){

    setwd(paste0(path, clu, "/"))
    file = list.files(path = ".", pattern = "_genotype-DEGs.csv")
    geno.deg.all[[clu]] <- read.csv(file)

}

int.deg.all <- list()
for (clu in cell.types){

    setwd(paste0(path, clu, "/"))
    file = list.files(path = ".", pattern = "_int-DEGs.csv")
    int.deg.all[[clu]] <- read.csv(file)

}

deg.all <- list()
for (clu in cell.types){

    thc = thc.deg.all[[clu]]
    thc$Type = rep("Treatment", dim(thc)[1])
    thc$Type[which(thc$padj >= 0.05)] = "NS"
    thc <- thc %>% distinct(Type, X, .keep_all = TRUE)
    geno = geno.deg.all[[clu]]
    geno$Type = rep("Genotype", dim(geno)[1])
    geno$Type[which(geno$padj >= 0.05)] = "NS"
    geno <- geno %>% distinct(Type, X, .keep_all = TRUE)
    int= int.deg.all[[clu]]
    int$Type = rep("Interaction", dim(int)[1])
    int$Type[which(int$padj >= 0.05)] = "NS"
    int <- int %>% distinct(Type, X, .keep_all = TRUE)
    deg.all[[clu]] <- rbind(thc, geno, int)

}

color_map <- c(
  Treatment   = "#E64B35",
  Genotype    = "#4DBBD5",
  Interaction = "#00A087",
  NS          = "grey75"
)

pics <- list()
for (clu in cell.types){

    df = deg.all[[clu]]
    col_custom <- setNames(color_map[df$Type], df$X)
    df$Type <- factor(df$Type,
                          levels = c("Treatment", "Genotype", "Interaction", "NS"))
    label_df <- df |>
      filter(Type != "NS") |>
      group_by(Type) |>
      slice_min(order_by = padj, n = 10) |>
      ungroup()

     p <- ggplot(df, aes(x = log2FC, y = -log10(padj), color = Type)) +
            geom_point(data = filter(df, Type == "NS"),
                        size = 1, alpha = 0.4) +
             geom_point(data = filter(df, Type != "NS"),
                        size = 1.5, alpha = 0.8) +
             geom_hline(yintercept = -log10(0.05),
                        linetype = "dashed", color = "grey40", linewidth = 0.4) +
             geom_vline(xintercept = c(-0.25, 0.25),
                        linetype = "dashed", color = "grey40", linewidth = 0.4) +
             geom_text_repel(
               data          = label_df,
               aes(label     = X),
               size          = 2.8,
               max.overlaps  = 20,
               box.padding   = 0.3,
               point.padding = 0.2,
               segment.color = "grey50",
               segment.size  = 0.3,
               show.legend   = FALSE
             ) +
             scale_color_manual(
               values = color_map,
               labels = c(
                 treatment   = "Treatment DEGs",
                 genotype    = "Genotype DEGs",
                 interaction = "Interaction DEGs",
                 NS          = "Not significant"
              )
             ) +
             labs(
               title = clu,
               x     = expression(log[2]~"Fold Change"),
               y     = expression(-log[10]~"(p"[adj]*")"),
               color = NULL
            ) +
             theme_classic(base_size = 13) +
             theme(
               legend.position  = "top",
               legend.key.size  = unit(0.4, "cm"),
               plot.title       = element_text(hjust = 0.5)
             ) +
             guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

    pics[[clu]] <- p

}

p <- Reduce( `+`, pics ) +
        patchwork::plot_layout(ncol = 5)
ggplot2::ggsave("20260611_volcano-DESeq2-broad.pdf", plot = p, width = 20, height = 10, units = "in")


## Overlap between single-cell DEGs and DESeq2 pseudobulk

load("./Merged/DEGs/20230803_DEGs/DEG.RData")
names(reln.markers)[which(names(reln.markers) == "GABA-1 (Cdh18+)")] <- "In-1"
names(reln.markers)[which(names(reln.markers) == "GABA-2 (Tiam2+)")] <- "NB-1"
names(reln.markers)[which(names(reln.markers) == "GABA-3 (St18+ Kit+)")] <- "In-2"
names(reln.markers)[which(names(reln.markers) == "GABA-4 (Cdh9+)")] <- "In-3"
names(reln.markers)[which(names(reln.markers) == "GABA-5 (Vip+ Cck+)")] <- "In-4"
names(reln.markers)[which(names(reln.markers) == "GABA-6 (Nos1+ Npy+)")] <- "In-5"
names(reln.markers)[which(names(reln.markers) == "GABA-7 (Calb2+)")] <- "NB-2"
names(reln.markers)[which(names(reln.markers) == "GABA-8 (Chat+)")] <- "In-6"
names(thc.markers)[which(names(thc.markers) == "GABA-1 (Cdh18+)")] <- "In-1"
names(thc.markers)[which(names(thc.markers) == "GABA-2 (Tiam2+)")] <- "NB-1"
names(thc.markers)[which(names(thc.markers) == "GABA-3 (St18+ Kit+)")] <- "In-2"
names(thc.markers)[which(names(thc.markers) == "GABA-4 (Cdh9+)")] <- "In-3"
names(thc.markers)[which(names(thc.markers) == "GABA-5 (Vip+ Cck+)")] <- "In-4"
names(thc.markers)[which(names(thc.markers) == "GABA-6 (Nos1+ Npy+)")] <- "In-5"
names(thc.markers)[which(names(thc.markers) == "GABA-7 (Calb2+)")] <- "NB-2"
names(thc.markers)[which(names(thc.markers) == "GABA-8 (Chat+)")] <- "In-6"

names(sig.reln.markers)[which(names(sig.reln.markers) == "GABA-1 (Cdh18+)")] <- "In-1"
names(sig.reln.markers)[which(names(sig.reln.markers) == "GABA-2 (Tiam2+)")] <- "NB-1"
names(sig.reln.markers)[which(names(sig.reln.markers) == "GABA-3 (St18+ Kit+)")] <- "In-2"
names(sig.reln.markers)[which(names(sig.reln.markers) == "GABA-4 (Cdh9+)")] <- "In-3"
names(sig.reln.markers)[which(names(sig.reln.markers) == "GABA-5 (Vip+ Cck+)")] <- "In-4"
names(sig.reln.markers)[which(names(sig.reln.markers) == "GABA-6 (Nos1+ Npy+)")] <- "In-5"
names(sig.reln.markers)[which(names(sig.reln.markers) == "GABA-7 (Calb2+)")] <- "NB-2"
names(sig.reln.markers)[which(names(sig.reln.markers) == "GABA-8 (Chat+)")] <- "In-6"
names(sig.thc.markers)[which(names(sig.thc.markers) == "GABA-1 (Cdh18+)")] <- "In-1"
names(sig.thc.markers)[which(names(sig.thc.markers) == "GABA-2 (Tiam2+)")] <- "NB-1"
names(sig.thc.markers)[which(names(sig.thc.markers) == "GABA-3 (St18+ Kit+)")] <- "In-2"
names(sig.thc.markers)[which(names(sig.thc.markers) == "GABA-4 (Cdh9+)")] <- "In-3"
names(sig.thc.markers)[which(names(sig.thc.markers) == "GABA-5 (Vip+ Cck+)")] <- "In-4"
names(sig.thc.markers)[which(names(sig.thc.markers) == "GABA-6 (Nos1+ Npy+)")] <- "In-5"
names(sig.thc.markers)[which(names(sig.thc.markers) == "GABA-7 (Calb2+)")] <- "NB-2"
names(sig.thc.markers)[which(names(sig.thc.markers) == "GABA-8 (Chat+)")] <- "In-6"

hg.test <- function(genes, marker.set, total.genes){
	      ngenes <- length(genes) # interested gene list
	      nmarker <- length(marker.set) # gene list being tested against for overlap
	      ntotal = length(unique(c(total.genes, marker.set)))
	      shared <- length(intersect(genes, marker.set))
	      p.value <- phyper(shared-1, ngenes, ntotal-ngenes, nmarker, lower.tail=FALSE)
	      return(p.value)
}

total.genes = data@assays[["SCT"]]@counts@Dimnames[[1]]

markers.list = list(sig.thc.markers, sig.reln.markers)
names(markers.list) = c("THC", "Reln")

result = list()
for (j in 1:length(markers.list)){

    n = markers.list[[j]]
    deg.overlap = data.frame(
                    Seurat.DEG.all = rep(0,length(names(n))),
                    Seurat.DEG.all.n = rep(0,length(names(n))),
                    DESeq2.DEG.all = rep(0, length(names(n))),
                    DESeq2.DEG.all.n = rep(0, length(names(n))),
                    Overlap.all = rep(0, length(names(n))),
                    Overlap.all.n = rep(0, length(names(n))),
                    HG.all = rep(0, length(names(n))),

                    Seurat.DEG.up = rep(0,length(names(n))),
                    Seurat.DEG.up.n = rep(0,length(names(n))),
                    DESeq2.DEG.up = rep(0, length(names(n))),
                    DESeq2.DEG.up.n = rep(0, length(names(n))),
                    Overlap.up = rep(0, length(names(n))),
                    Overlap.up.n = rep(0, length(names(n))),
                    HG.up = rep(0, length(names(n))),

                    Seurat.DEG.down = rep(0,length(names(n))),
                    Seurat.DEG.down.n = rep(0,length(names(n))),
                    DESeq2.DEG.down = rep(0, length(names(n))),
                    DESeq2.DEG.down.n = rep(0, length(names(n))),
                    Overlap.down = rep(0, length(names(n))),
                    Overlap.down.n = rep(0, length(names(n))),
                    HG.down = rep(0, length(names(n)))
                    )
    rownames(deg.overlap) = names(n)

    for(i in 1:length(names(n))){
        temp = n[[i]]
        sc.deg.up = rownames(temp)[which(temp$avg_log2FC > 0)]
        sc.deg.down = rownames(temp)[which(temp$avg_log2FC < 0)]

        deg.overlap$Seurat.DEG.all[i] = paste(rownames(temp), collapse = ", ")
        deg.overlap$Seurat.DEG.all.n[i] = dim(temp)[1]
        deg.overlap$Seurat.DEG.up[i] = paste(sc.deg.up, collapse = ", ")
        deg.overlap$Seurat.DEG.up.n[i] = length(sc.deg.up)
        deg.overlap$Seurat.DEG.down[i] = paste(sc.deg.down, collapse = ", ")
        deg.overlap$Seurat.DEG.down.n[i] = length(sc.deg.down)

        if(j == 1){
            temp.ps = thc.deg[[names(n)[i]]]
        }else{
            temp.ps = geno.deg[[names(n)[i]]]
        }
        ps.deg.up = temp.ps$X[which(temp.ps$log2FC > 0)]
        ps.deg.down = temp.ps$X[which(temp.ps$log2FC < 0)]

        deg.overlap$DESeq2.DEG.all[i] = paste(temp.ps$X, collapse = ", ")
        deg.overlap$DESeq2.DEG.all.n[i] = length(temp.ps$X)
        deg.overlap$DESeq2.DEG.up[i] = paste(ps.deg.up, collapse = ", ")
        deg.overlap$DESeq2.DEG.up.n[i] = length(ps.deg.up)
        deg.overlap$DESeq2.DEG.down[i] = paste(ps.deg.down, collapse = ", ")
        deg.overlap$DESeq2.DEG.down.n[i] = length(ps.deg.down)

        overlap.deg = intersect(rownames(temp), temp.ps$X)
        deg.overlap$Overlap.all[i] = paste(overlap.deg, collapse = ", ")
        deg.overlap$Overlap.all.n[i] = length(overlap.deg)
        deg.overlap$HG.all[i] = hg.test(rownames(temp), temp.ps$X, total.genes)

        overlap.up = intersect(sc.deg.up, ps.deg.up)
        deg.overlap$Overlap.up[i] = paste(overlap.up, collapse = ", ")
        deg.overlap$Overlap.up.n[i] = length(overlap.up)
        deg.overlap$HG.up[i] = hg.test(sc.deg.up, ps.deg.up, total.genes)

        overlap.down = intersect(sc.deg.down, ps.deg.down)
        deg.overlap$Overlap.down[i] = paste(overlap.down, collapse = ", ")
        deg.overlap$Overlap.down.n[i] = length(overlap.down)
        deg.overlap$HG.down[i] = hg.test(sc.deg.down, ps.deg.down, total.genes)

    }
    result[[names(markers.list)[j]]] <- deg.overlap

}
WriteXLS(result, "./Merged/DEGs/20260504_DESeq2/Covariate_depth/Seurat-DESeq2-DEG-overlap.xlsx", names(result), row.names = TRUE)
