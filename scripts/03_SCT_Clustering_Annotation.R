setwd("~/Desktop/HRM-THC-snRNA-seq/Merged")
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
