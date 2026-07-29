setwd("./Merged")
library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)
library(patchwork)
library(here)
set.seed(12)
# parallelization
library(future)
plan("multicore", workers = 8)
plan()
options(future.globals.maxSize = 100 * 1024 ^ 3) # for 100 Gb RAM
# run sctransform
# regressing out nCount_RNA, percent.mito, percent.rpl, percent.rps
data <- SCTransform(data, vars.to.regress = c("percent.mt", "nCount_RNA", "percent.rpl", "percent.rps"), verbose = TRUE)
save.image("20230104_68917_SCT-reg-full.RData")
# clustering
DefaultAssay(data) <- "SCT"
data <- RunPCA(data, verbose = TRUE)
ElbowPlot(data, ndims = 50) # to determine number of dimensions for clustering
data <- FindNeighbors(data, dims = 1:50, verbose = TRUE)
data <- FindClusters(data, verbose = TRUE, resolution = 0.7)
data <- RunUMAP(data, dims = 1:50, verbose = TRUE)
# cell type annotation
p1 <- DimPlot(data, reduction = "umap", label = TRUE) + ggtitle("snRNA Seurat Clustering")
celltype.markers <- c("Snap25", "Syt1", "Camk2a", # neuron
                      "Slc17a6", "Slc17a7", # excitatory
                      "Resp18", "Gad1", "Gad2", # inhibitory
                      "Drd1", "Pdyn", "Tac1", # D1-MSN
                      "Drd2","Adora2a", "Penk", # D2-MSN
                      "Grm8", "Drd3", "Ppp1r1b", "Foxp2", "Bcl11b", # DA
                      "Mbp", # Ol
                      "Pdgfra", # OPC
                      "Gja1", "Aldoc", "Aqp4", # Ast
                      "Top2a", # NB neural stem cells and neuroblasts
                      "C1qa", "C1qc", "Cd74", "Arhgap15", # Mg
                      "Cd163", # Macrophage
              "Vtn", # Mural
              "Cldn5", # Endo
              "Ccdc153") # Ependymal
p2 <- DotPlot(data, features = celltype.markers) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
CombinePlots(plots = list(p1, p2))
data[["orig.clusters"]] <- Idents(object = data) # stash cluster idents prior to annotation
# double check some markers
markers <- list()
for (clu in levels(Idents(data))){
    markers[[clu]] <- FindMarkers(object = data, ident.1 = clu)
}
# detailed annotation
new.cluster.ids <- c("D1.MSN-1","D2.MSN-1","Ol-1", "Ast-1", "D1.MSN-2",
                     "D2.MSN-2","D1.MSN-3","Glut-1","D1.MSN-4","D1.MSN-5",
                     "GABA-1","D2.MSN-3","GABA-2","OPC","Glut-2",
                     "Glut-3","D3.MSN", "GABA-3", "GABA-4", "Ast-2",
                     "Mg", "Ol-2", "Glut-4", "Ast-3", "GABA-5",
                     "GABA-6", "26", "Glut-5", "GABA-7", "NB", "Mu-Endo", "GABA-8",
                     "Glut-6", "Ol-3")
names(new.cluster.ids) <- levels(data)
data <- RenameIdents(data, new.cluster.ids)
levels(data) <- c("D1.MSN-1","D1.MSN-2", "D1.MSN-3","D1.MSN-4","D1.MSN-5","D3.MSN",
               "D2.MSN-1","D2.MSN-2","D2.MSN-3",
               "GABA-1","GABA-2","GABA-3","GABA-4","GABA-5","GABA-6","GABA-7","GABA-8",
               "Glut-1","Glut-2","Glut-3","Glut-4","Glut-5","Glut-6",
               "Ol-1","Ol-2","Ol-3","OPC",
               "Ast-1","Ast-2","Ast-3","NB", "Mg", "Mu-Endo","26")
data@meta.data$detailed.celltype <- data@active.ident
names(markers) <- new.cluster.ids
# remove cluster 26 due to high level of mito genes in markers
# remove cluster 33 due to mixed markers
FeaturePlot(data, features = "percent.mt") + ggtitle("Mito.genes distribution")
VlnPlot(data, features = c("percent.mt", "percent.rpl", "percent.rps"), ncol = 1, pt.size=0)
data <- subset(x = data, idents = c("26", "Ol-3"), invert = TRUE)
levels(data) <- c("D1.MSN-1","D1.MSN-2", "D1.MSN-3","D1.MSN-4","D1.MSN-5","D3.MSN",
               "D2.MSN-1","D2.MSN-2","D2.MSN-3",
               "GABA-1","GABA-2","GABA-3","GABA-4","GABA-5","GABA-6","GABA-7","GABA-8",
               "Glut-1","Glut-2","Glut-3","Glut-4","Glut-5","Glut-6",
               "Ol-1","Ol-2","OPC",
               "Ast-1","Ast-2","Ast-3","NB", "Mg", "Mu-Endo")
data@meta.data$detailed.celltype <- data@active.ident
# coarse annotation
new.cluster.ids <- c("D1.MSN","D1.MSN", "D1.MSN","D1.MSN","D1.MSN","D3.MSN",
               "D2.MSN","D2.MSN","D2.MSN",
               "GABA","GABA","GABA","GABA","GABA","GABA","GABA","GABA",
               "Glut","Glut","Glut","Glut","Glut","Glut",
               "Ol","Ol","OPC",
               "Ast","Ast","Ast","NB", "Mg", "Mu-Endo")
names(new.cluster.ids) <- levels(data)
data <- RenameIdents(data, new.cluster.ids)
levels(data) <- c("D1.MSN","D3.MSN",
               "D2.MSN",
               "GABA",
               "Glut",
               "Ol","OPC",
               "Ast","NB", "Mg", "Mu-Endo")
data@meta.data$coarse.celltype <- data@active.ident
DimPlot(data, reduction = "umap", label = TRUE) + ggtitle("Coarse cell types")
# Remove Glut population based on previous literature
Idents(object = data) <- data@meta.data[["coarse.celltype"]]
data <- subset(x = data, idents = "Glut", invert = TRUE) # n = 61911
save.image("20230107_SCT_noGlut-61911.RData")
## renaming
new.cluster.ids <- c("D1.MSN-1", "D1.MSN-2", "D1.MSN-3", "D1.MSN-4", "D1.MSN-5",
               "D2.MSN-1", "D2.MSN-2", "D2.MSN-3", "D3.MSN",
               "In-1", "NB-1", "In-2", "In-3", "In-4", "In-5", "NB-2", "In-6", "NPC",
               "Ol-1", "Ol-2", "OPC",
               "Ast-1", "Ast-2", "Ast-3", "Mg", "Mu-Endo")
names(new.cluster.ids) <- levels(data)
data <- RenameIdents(data, new.cluster.ids)
levels(data) <- c("D1.MSN-1", "D1.MSN-2", "D1.MSN-3", "D1.MSN-4", "D1.MSN-5",
               "D2.MSN-1", "D2.MSN-2", "D2.MSN-3", "D3.MSN",
               "In-1", "In-2", "In-3", "In-4", "In-5", "In-6",
               "NB-1", "NB-2", "NPC",
               "Ol-1", "Ol-2", "OPC",
               "Ast-1", "Ast-2", "Ast-3", "Mg", "Mu-Endo")
data@meta.data$detailed.celltype <- data@active.ident
save.image(file = "20250324.RData")

# ============================================================
# Visualization: cell-type marker violin plots, UMAP, cell counts
# ============================================================

setwd("./Merged")

library(cowplot)

color <- read.csv("cell-color-palette.csv")

celltype.markers <- c("Gad1", "Bcl11b", "Ppp1r1b", "Foxp2", "Drd1", "Drd2", "Drd3", "Htr4", "Ebf1", "Scube1", "Stk32a", "Grm8", "Chst9", "Oprm1", "Kit", "Sst", "Vip", "Cck", "Nos1", "Npy", "Calb2", "Chat", "Reln", "Top2a", "Gja1", "Mbp", "Pdgfra", "C1qa", "Cldn5", "Vtn")

levels(data) <- rev(c("D1.MSN-1","D1.MSN-2", "D1.MSN-3","D1.MSN-4","D1.MSN-5",               "D2.MSN-1","D2.MSN-2","D2.MSN-3", "D3.MSN",
               "GABA-1 (Cdh18+)","GABA-2 (Igfbpl1+)","GABA-3 (St18+ Kit+)","GABA-4 (Cdh9+)","GABA-5 (Vip+ Cck+)","GABA-6 (Nos1+ Npy+)","GABA-7 (Calb2+)","GABA-8 (Chat+)", "NPC",               "Ol-1","Ol-2","OPC",               "Ast-1","Ast-2","Ast-3", "Mg", "Mu-Endo"))

VlnPlot(data, celltype.markers, stack = TRUE) +
        theme(legend.position = "none") + labs(y = "")

# fetch MSN subtype markers
# NOTE: this snapshot labels the GABA-2 cluster by its Tiam2+ marker, unlike the
# Igfbpl1+ label used in the previous VlnPlot block above (from `data` as it
# stood after the earlier renaming steps) -- the levels<- call below only
# resolves if `data` actually comes from this different snapshot, hence the load.
load("./Merged/20230910.RData")

levels(data) <- rev(c("D1.MSN-1","D1.MSN-2", "D1.MSN-3","D1.MSN-4","D1.MSN-5",               "D2.MSN-1","D2.MSN-2","D2.MSN-3", "D3.MSN",
               "GABA-1 (Cdh18+)","GABA-2 (Tiam2+)","GABA-3 (St18+ Kit+)","GABA-4 (Cdh9+)","GABA-5 (Vip+ Cck+)","GABA-6 (Nos1+ Npy+)","GABA-7 (Calb2+)","GABA-8 (Chat+)", "NB",               "Ol-1","Ol-2","OPC",               "Ast-1","Ast-2","Ast-3", "Mg", "Mu-Endo"))

celltype.markers <- c("Gad1", "Bcl11b", "Ppp1r1b", "Foxp2", "Drd1", "Drd2", "Drd3", "Ebf1", "Onecut2", "Sema5b", "Greb1l", "Grm8", "Chst9", "Oprm1", "Slit2", "Rgs6", "Htr7", "Kit", "Sst", "Vip", "Cck", "Nos1", "Npy", "Calb2", "Chat", "Reln", "Top2a", "Mbp", "Pdgfra", "Gja1", "C1qa", "Cldn5", "Vtn")

VlnPlot(data, celltype.markers, stack = TRUE) +
        theme(legend.position = "none") + labs(y = "")

#### remove cell markers from idents

DefaultAssay(data) <- "SCT"
Idents(object = data) <- data@meta.data[["detailed.celltype"]]

new.cluster.ids <- c("D1.MSN-1","D1.MSN-2", "D1.MSN-3","D1.MSN-4","D1.MSN-5", "D3.MSN",               "D2.MSN-1","D2.MSN-2","D2.MSN-3",
               "GABA-1","GABA-2","GABA-3","GABA-4","GABA-5","GABA-6","GABA-7","GABA-8",                "Ol-1","Ol-2","OPC",               "Ast-1","Ast-2","Ast-3", "NPC", "Mg", "Mu-Endo")

names(new.cluster.ids) <- levels(data)
data <- RenameIdents(data, new.cluster.ids)

levels(data) <- c("D1.MSN-1","D1.MSN-2", "D1.MSN-3","D1.MSN-4","D1.MSN-5",               "D2.MSN-1","D2.MSN-2","D2.MSN-3", "D3.MSN",
               "GABA-1","GABA-2","GABA-3","GABA-4","GABA-5","GABA-6","GABA-7","GABA-8", "NPC",               "Ol-1","Ol-2","OPC",               "Ast-1","Ast-2","Ast-3", "Mg", "Mu-Endo")

data@meta.data$detailed.celltype <- data@active.ident
save.image(file = "20250130.RData")

# new dimplot
DimPlot(data, reduction = "umap", label = TRUE, repel = TRUE, cols = color$Color)

# vln
celltype.markers <- c("Gad1", "Ppp1r1b", "Foxp2", "Drd1", "Drd2", "Drd3", "Tac1", "Adora2a", "Ebf1", "Onecut2", "Sema5b", "Greb1l", "Grm8", "Htr7", "Resp18", "Kit", "Sst", "Crhbp", "Vip", "Cck", "Nos1", "Npy", "Calb2", "Chat", "Top2a", "Mbp", "Pdgfra", "Gja1", "C1qa", "Cldn5", "Vtn")

VlnPlot(data, celltype.markers, stack = TRUE) +
        theme(legend.position = "none") + labs(y = "")

# vertical
levels(data) <- c("D1.MSN-1","D1.MSN-2", "D1.MSN-3","D1.MSN-4","D1.MSN-5",               "D2.MSN-1","D2.MSN-2","D2.MSN-3", "D3.MSN",
               "GABA-1","GABA-2","GABA-3","GABA-4","GABA-5","GABA-6","GABA-7","GABA-8", "NPC",               "Ol-1","Ol-2","OPC",               "Ast-1","Ast-2","Ast-3", "Mg", "Mu-Endo")

VlnPlot(data, celltype.markers, stack = TRUE, flip = TRUE) +
        theme(legend.position = "none") + labs(x = "")

celltype.markers <- c("Gad1", "Reln", "Cnr1")
VlnPlot(data, celltype.markers, stack = TRUE, flip = TRUE) +
        theme(legend.position = "none") + labs(x = "")

celltype.markers <- c("Gad1", "Reln", "Cnr1", "Lrp8")
VlnPlot(data, celltype.markers, stack = TRUE, flip = TRUE) +
        theme(legend.position = "none") + labs(x = "")

# feature plot on the original UMAP

load("./Merged/Archive/20230104_68917_SCT-reg-full.RData")

DefaultAssay(data) <- "SCT"
data <- RunPCA(data, verbose = TRUE)
data <- FindNeighbors(data, dims = 1:50, verbose = TRUE)
data <- FindClusters(data, verbose = TRUE, resolution = 0.7)
data <- RunUMAP(data, dims = 1:50, verbose = TRUE)

p1 <- DimPlot(data, reduction = "umap", label = TRUE, repel = TRUE)
p2 <- FeaturePlot(data, features = "Slc17a7")
p1 + p2

# Cell count barplot

df <- data.frame(table(Idents(data)))
colnames(df) <- c("Cell type", "Counts")
df$Color <- color$Color

df$`Cell type` <- factor(df$`Cell type`, levels = rev(df$`Cell type`))

ggplot(df, aes(fill = `Cell type`, x = Counts, y =`Cell type`)) +    geom_bar(stat = "identity", show.legend = FALSE) + theme_classic() +    labs(y = "", x = "Number of nuclei") +    theme(text = element_text(size = 12)) +    geom_text(aes(label = Counts), vjust = -0.2) +     scale_fill_manual(values = rev(df$Color)) +
    theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))

# adjusted feature plot

celltype.markers <- c("Gad1", # inhibitory
                      "Foxp2",
                      "Drd1",  # D1-MSN
                      "Drd2", # D2-MSN
                      "Drd3",  # DA
                      "Resp18", # interneuron
                      "Top2a", # NB neural stem cells and neuroblasts
                      "Mbp", # Ol
                      "Pdgfra", # OPC
                      "Gja1", # Ast
                      "C1qa", # Mg
		      "Vtn") # Mural

FeaturePlot(data, features = celltype.markers, ncol = 3)
