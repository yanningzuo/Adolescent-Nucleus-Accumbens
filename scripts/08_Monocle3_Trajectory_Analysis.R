library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(patchwork)
library(dplyr)

set.seed(12)

setwd("./Merged")

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

# ============================================================
# Visualization: trajectory-branch DEG enrichment + cell-composition
# pie charts (originally Scripts/20231015_traj1-vis.txt)
# ============================================================

deg <- read.csv("traj1-pseudotime-deg.csv")
deg <- deg[which(deg$q_value < 0.05), ]

library(enrichR)
dbs <- c("KEGG_2021_Human", "GO_Biological_Process_2021", "BioCarta_2016", "Reactome_2022")

genes <- deg$X[1:100]
go.genes <- enrichr(genes, dbs)
list.go.genes = do.call("rbind",go.genes)
list.go.genes$padj <- p.adjust(list.go.genes$P.value, method = "BH")
list.go.genes <- list.go.genes[order(list.go.genes$padj), ]

## scCODA Prep

DefaultAssay(traj1.neurons) <- "SCT"
Idents(object = traj1.neurons) <- traj1.neurons@meta.data[["orig.ident"]]

cell_counts = ''

for (sample in unique(traj1.neurons@meta.data[["orig.ident"]])){
    temp <- subset(x = traj1.neurons, subset = orig.ident == sample)
    Idents(object = temp) <- temp@meta.data[["subclusters"]]
    df <- data.frame(table(temp@meta.data[["subclusters"]]))
    cell_counts = as.data.frame(cbind(cell_counts, df[, 2]))
    rownames(cell_counts) = df[, 1]; colnames(cell_counts)[dim(cell_counts)[2]] = sample
}

cell_counts <- cell_counts[, -1]
cell_counts <- t(cell_counts)
write.csv(cell_counts, file = "traj1-scCODA_cell_counts.csv", row.names = T)

# pie charts of treatment percentage per cluster

Idents(object = traj1.neurons) <- traj1.neurons@meta.data[["Treatment"]]
df <- data.frame(table(traj1.neurons@active.ident))
colnames(df) <- c("Group", "value")

# df$Group <- factor(df$Group, levels = c("THC", "Veh"))

# Compute the position of labels
df <- df %>%
  arrange(desc(Group)) %>%
  mutate(prop = value / sum(df$value) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

# Basic piechart
p1 <- ggplot(df, aes(x="", y=prop, fill=Group)) +
  	geom_bar(stat="identity", width=1, color="white") +
  	coord_polar("y", start=0) +
  	theme_void() +
  	#geom_text(aes(y = ypos, label = paste0(round(prop, 2), "%")), color = "white", size=6)+
  	ggtitle("All cells")

pics <- list()
pics[["all"]] <- p1
for (clu in levels(traj1.neurons)){
    temp <- subset(x = traj1.neurons, subset = subclusters == clu)
    Idents(object = temp) <- temp@meta.data[["Treatment"]]
    df <- data.frame(table(temp@active.ident))
    colnames(df) <- c("Group", "value")
    df <- df %>%
      arrange(desc(Group)) %>%
      mutate(prop = value / sum(df$value) *100) %>%
      mutate(ypos = cumsum(prop)- 0.5*prop )
    p <- ggplot(df, aes(x="", y=prop, fill=Group)) +
      	    geom_bar(stat="identity", width=1, color="white") +
  	    coord_polar("y", start=0) +
  	    theme_void() +
  	    ggtitle(clu)
    pics[[clu]] <- p
}

Reduce( `+`, pics ) +
      patchwork::plot_layout(ncol = 6)

# pie charts of genotype percentage per cluster

Idents(object = traj1.neurons) <- traj1.neurons@meta.data[["Genotype"]]
df <- data.frame(table(traj1.neurons@active.ident))
colnames(df) <- c("Group", "value")

# Compute the position of labels
df <- df %>%
  arrange(desc(Group)) %>%
  mutate(prop = value / sum(df$value) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

# Basic piechart
p1 <- ggplot(df, aes(x="", y=prop, fill=Group)) +
  	geom_bar(stat="identity", width=1, color="white") +
  	coord_polar("y", start=0) +
  	theme_void() +
  	#geom_text(aes(y = ypos, label = paste0(round(prop, 2), "%")), color = "white", size=6)+
  	ggtitle("All cells")

pics <- list()
pics[["all"]] <- p1
for (clu in unique(traj1.neurons@meta.data[["subclusters"]])){
    temp <- subset(x = traj1.neurons, subset = subclusters == clu)
    Idents(object = temp) <- temp@meta.data[["Genotype"]]
    df <- data.frame(table(temp@meta.data[["Genotype"]]))
    colnames(df) <- c("Group", "value")
    df <- df %>%
      arrange(desc(Group)) %>%
      mutate(prop = value / sum(df$value) *100) %>%
      mutate(ypos = cumsum(prop)- 0.5*prop )
    p <- ggplot(df, aes(x="", y=prop, fill=Group)) +
      	    geom_bar(stat="identity", width=1, color="white") +
  	    coord_polar("y", start=0) +
  	    theme_void() +
  	    ggtitle(clu)
    pics[[clu]] <- p
}

Reduce( `+`, pics ) +
      patchwork::plot_layout(ncol = 6)
