library(SCENIC)
library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)
library(patchwork)

set.seed(12)

# parallelization
library(future)
plan("multicore", workers = 12)
plan()
options(future.globals.maxSize = 200 * 1024 ^ 3) # for 100 Gb RAM

setwd("./Merged")
load("20230430_61911.RData")

# Data prep
DefaultAssay(data) <- "RNA"
Idents(object = data) <- data@meta.data[["coarse.celltype"]]
group.cells <- WhichCells(object = data, idents = "GABA")
exprMat <- as.matrix(x = data@assays[["RNA"]]@data[, group.cells])
Idents(object = data) <- data@meta.data[["detailed.celltype"]]
cellInfo <- data.frame(seuratCluster = Idents(data))

# Initialize settings
setwd(".")
org <- "mgi"
dbDir <- "cisTarget_databases"
myDatasetTitle <- "GABA-all"
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
lapply(dbs,function(x) file.exists(file.path(dbDir, x)))

saveRDS(cellInfo, file="int/cellInfo.Rds")

motifAnnotations_mgi <- motifAnnotations
scenicOptions <- initializeScenic(org = "mgi", dbDir = "cisTarget_databases", datasetTitle=myDatasetTitle, nCores = 12)
motifAnnotations_mgi <- motifAnnotations
scenicOptions@settings$seed <- 12

saveRDS(scenicOptions, file="int/scenicOptions.Rds")

# Co-expression network
genesKept <- geneFiltering(exprMat, scenicOptions, minCountsPerGene=3*.01*ncol(exprMat), minSamples=ncol(exprMat)*.01)
exprMat_filtered <- exprMat[genesKept, ]

runCorrelation(exprMat_filtered, scenicOptions)

exprMat_filtered_log <- log2(exprMat_filtered+1)
runGenie3(exprMat_filtered_log, scenicOptions)

# Build and score the GRN
setwd(".")
exprMat_log <- log2(exprMat+1)

scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 12
scenicOptions@settings$seed <- 12

scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, minGenes = 20)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)

export2loom(scenicOptions, exprMat)
save.image(file = paste0("scenic-", myDatasetTitle, ".RData"))

library(SCopeLoomR)

setwd("./Archive/Scenic/scenic-gaba-all-20231210")
loom <- open_loom("./Archive/Scenic/scenic-gaba-all-20231210/output/scenic.loom")

# Read information from loom file:
regulons_incidMat <- get_regulons(loom)
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonsAUC <- get_regulons_AUC(loom)
regulonsAucThresholds <- get_regulon_thresholds(loom)
regulonsAUC <- regulonsAUC[onlyNonDuplicatedExtended(rownames(regulonsAUC)),]
embeddings <- get_embeddings(loom)

load("./Archive/Scenic/scenic-gaba-all-20231210/scenic-GABA-all.RData")
Idents(object = data) <- data@meta.data[["detailed.celltype"]]
cellInfo <- data.frame(seuratCluster=Idents(data))
colnames(cellInfo) <- 'CellType'
cellInfo$barcode = rownames(cellInfo)
gaba = c("In-1", "In-2", "In-3", "In-4", "In-5", "In-6", "NB-1", "NB-2", "NPC")
cellInfo <- cellInfo[which(cellInfo$CellType %in% gaba), ]
cellInfo$barcode <- NULL

regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                     function(cells) rowMeans(getAUC(regulonsAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
regulonActivity_byCellType_Scaled <- regulonActivity_byCellType_Scaled[, 10:17]

library(ComplexHeatmap)
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity")

topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]

rss <- calcRSS(AUC=getAUC(regulonsAUC), cellAnnotation=cellInfo[colnames(regulonsAUC), "CellType"])
write.csv(rss, "rss.csv")
