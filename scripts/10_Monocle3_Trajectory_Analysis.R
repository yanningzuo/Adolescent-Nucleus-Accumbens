library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(patchwork)
library(dplyr)

set.seed(12)

setwd("~/Desktop/HRM-THC-snRNA-seq/Merged")

DefaultAssay(data) <- "SCT"
Idents(object = data) <- data@meta.data[["detailed.celltype"]]

## trajectory inference with interneunons, neuroblasts, and NPCs

nb.int <- subset(x = data, idents = c("In-1", "In-2", "In-3", "In-4", "In-5", "In-6", "NB-1", "NB-2", "NPC"))
cds <- as.cell_data_set(nb.int)
fData(cds)$gene_short_name <- rownames(fData(cds))

recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
recreate.partitions

cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions

list.cluster <- nb.int@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster

cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- nb.int@reductions$umap@cell.embeddings

cluster.before.traj <- plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F, 
           group_label_size = 5) + theme(legend.position = "right")
cluster.before.traj

cds <- learn_graph(cds, use_partition = F)
p1 <- plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F,
           label_branch_points = T, label_roots = T, label_leaves = F,
           group_label_size = 5)

cds <- order_cells(cds, reduction_method = "UMAP", root_cells = colnames(cds[, clusters(cds) == "NB"]))
p2 <- plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = T,
           label_branch_points = T, label_roots = F, label_leaves = F)

p1 + p2

cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))
data.pseudo$detailed.celltype <- factor(data.pseudo$detailed.celltype, levels = rev(c("In-1", "In-2", "In-3", "In-4", "In-5", "In-6", "NB-1", "NB-2", "NPC")))

p3 <- ggplot(data.pseudo, aes(monocle3_pseudotime, detailed.celltype, fill = detailed.celltype)) + geom_boxplot()
(p1 | p2) /
p3 

## trajectory inference with all neuronal cell types

data <- subset(x = data, idents = levels(data)[c(1:17, 24)])
cds <- as.cell_data_set(data)
cds <- cluster_cells(cds, resolution = 1e-5)
fData(cds)$gene_short_name <- rownames(fData(cds))

list.cluster <- data@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster

cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- data@reductions$umap@cell.embeddings
plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F, 
           group_label_size = 5) + theme(legend.position = "right")
plot_cells(cds, color_cells_by="partition", group_cells_by="cluster")

cds <- learn_graph(cds, use_partition = F)
p1 <- plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F,
           label_branch_points = T, label_roots = T, label_leaves = F,
           group_label_size = 5)

cds <- order_cells(cds, reduction_method = "UMAP", root_cells = colnames(cds[, clusters(cds) == "NB"]))
p2 <- plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = T,
           label_branch_points = T, label_roots = F, label_leaves = F)

p1 + p2

cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))
p3 <- ggplot(data.pseudo, aes(monocle3_pseudotime, detailed.celltype, fill = detailed.celltype)) + geom_boxplot()
(p1 | p2) /
p3
