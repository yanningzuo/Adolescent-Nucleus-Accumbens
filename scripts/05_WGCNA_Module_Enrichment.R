setwd("~/Desktop/HRM-THC-snRNA-seq/Merged")
load("~/Desktop/HRM-THC-snRNA-seq/Merged/20230123_SCT_noGlut-61911.RData")

library(Seurat)
set.seed(12)

# parallelization
library(future)
plan("multicore", workers = 8)
plan()

# prepare inputs for WGCNA

DefaultAssay(data) <- "RNA"
Idents(object = data) <- data@meta.data[["detailed.celltype"]]
data <- NormalizeData(object = data, normalization.method = "LogNormalize")

outpath = "~/Desktop/HRM-THC-snRNA-seq/WGCNA-RNA/Inputs/"

for(group in levels(data@active.ident)) {
    group.cells <- WhichCells(object = data, idents = group)
    data_to_write_out <- as.data.frame(as.matrix(x = data@assays[["RNA"]]@data[, group.cells]))
    write.csv(x = data_to_write_out, row.names = TRUE, file = paste0(outpath, group, "_RNA-log1p.csv"))
}

# scWGCNA

library(WGCNA)
library(flashClust)
options(stringsAsFactors = FALSE);
allowWGCNAThreads() 

path = c("/home/y3zuo/gpfs/HRM-THC/WGCNA/")
name <- c("HRM-THC-wgcna")
source <- "NAcc"
setwd(file.path(path, name))
input_paths = c("/home/y3zuo/gpfs/HRM-THC/WGCNA/")
files = list.files(input_paths)

for (f in files){
    filename <- gsub('_.*$', '', f)
    dir.create(file.path(path, source, filename))
    outpath = file.path(paste0(path, source, "/" , filename, "/"))

    # extraction of data and thresholding
    data <- read.csv(paste0(input_paths, f), header = TRUE, row.names = 1)
    dseurat <- CreateSeuratObject(counts = data)
    low.limit <- dseurat@meta.data[["nCount_RNA"]][match(round(length(dseurat@meta.data[["nCount_RNA"]])*0.4), rank(-dseurat@meta.data[["nCount_RNA"]], ties.method = c("first")))]
    high.limit <- dseurat@meta.data[["nCount_RNA"]][match(round(length(dseurat@meta.data[["nCount_RNA"]])*0.05), rank(-dseurat@meta.data[["nCount_RNA"]], ties.method = c("first")))]
    dseurat <- subset(dseurat, subset = nCount_RNA > low.limit  & nCount_RNA < high.limit)
    data <- as.data.frame(as.matrix(x = dseurat@assays[["RNA"]]@data))

    # WGCNA preparation
    sd = apply(data, 1, sd)
    data = data[sd > 0, ]
    data=t(data)
    gsg=goodSamplesGenes(data, verbose=3)
    gsg$allOK

    # sft calculation
    powers = c(c(1:11), seq(from = 12, to=30, by=2))
    sft = pickSoftThreshold(data, networkType = "signed",corFnc = "bicor",verbose = 5,powerVector = powers)
    if (length(which(sft$fitIndices[, 2] >= 0.8)) == 0){
        softPower <- sft$fitIndices[, 1][which.max(sft$fitIndices[, 2])]
    } else {
        softPower <- sft$fitIndices[, 1][min(which(sft$fitIndices[, 2] >= 0.8))] # was this - 1
    }
    write.csv(sft, file = paste0(outpath, filename, "_sft.csv"))

    # plot sft
    pdf(paste0(outpath, filename, "_", softPower, "_sft.pdf"))
    par(mfrow = c(1,2));
    cex1 = 0.9;
    plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type="n", main = paste("Scale independence"))
    text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], labels=powers, cex=cex1, col="red")
    abline(h = 0.90, col = "red")
    plot(sft$fitIndices[, 1], sft$fitIndices[, 5], xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n", main = paste("Mean connectivity"))
    text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")
    dev.off()

    # calculating modules
    adj = adjacency(data, type = "signed", power = softPower)
    TOM = TOMsimilarity(adj, TOMType= "signed", verbose = 0)
    dissTOM = 1 - TOM
    geneTree = flashClust(as.dist(dissTOM), method = "average")
    minModuleSize = 30
    ds = 2
    cutHeight = 0.9999
    pamStage = F
    dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, method = "hybrid", deepSplit = ds, pamRespectsDendro = T, pamStage = pamStage,
                            minClusterSize = minModuleSize, cutHeight = cutHeight)
    dynamicColors = labels2colors(dynamicMods)

    pdf(paste0(outpath, filename, "_module.pdf"))
    plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = filename)
    dev.off()

    # clusters
    MEList = moduleEigengenes(data, colors = dynamicColors, softPower = softPower)
    MEs = MEList$eigengenes
    MEDiss = 1 - cor(MEs)
    METree = hclust(as.dist(MEDiss), method = "average")

    pdf(paste0(outpath, filename, "_cluster.pdf"))
    plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
    MEDissThres = 0.15
    abline(h = MEDissThres, col = "red")
    dev.off()

    # merge modules
    merge = mergeCloseModules(data, dynamicColors, cutHeight = MEDissThres, verbose = 3)
    mergedColors = merge$colors

    pdf(paste0(outpath, filename, "_merged.pdf"))
    plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
    dev.off()

    moduleColors = mergedColors
    MEs = merge$newMEs
    colorOrder = c("grey", standardColors(50))
    moduleLabels = match(moduleColors, colorOrder) - 1
    modNames = substring(names(MEs), 3)
    totalmodnum <- length(unique(mergedColors))

    geneModulebicor = corAndPvalue(data, MEs, use = "p")
    geneModuleMembership = as.data.frame(geneModulebicor$cor)
    MMPvalue = as.data.frame(geneModulebicor$p)
    names(geneModuleMembership) = paste("MM", modNames, sep = "")
    names(MMPvalue) = paste("p.MM", modNames, sep = "")

    # Output
    geneInfo0 = data.frame(GENE = colnames(data), moduleColor = moduleColors)
    for (mod in 1:ncol(geneModuleMembership)){
        oldNames = names(geneInfo0)
        geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, mod], MMPvalue[, mod])
        names(geneInfo0) = c(oldNames, paste("MM.", modNames[mod], sep = ""), paste("p.MM.", modNames[mod], sep = ""))
    }
    geneOrder = order(geneInfo0$moduleColor)
    geneInfo = geneInfo0[geneOrder, ]
    geneInfo$MODULE = 0
    count = 1
    for (i in 1:(dim(geneInfo)[1]-1)){
        if (geneInfo$moduleColor[i] == geneInfo$moduleColor[i+1]){
            geneInfo$MODULE[i] = count
        } else {
            geneInfo$MODULE[i] = count
            count = count + 1
        }
        geneInfo$MODULE[dim(geneInfo)[1]] = dim(MEs)[2]
    }

    output <- geneInfo[, c("GENE","MODULE","moduleColor")]
    output$moduleColor <- gsub("darkgrey", "darkgray", output$moduleColor)
    output$moduleColor <- gsub("grey60", "gray60", output$moduleColor)
    output <- output[!grepl("grey", output$moduleColor), ]
    write.csv(geneInfo, file = paste0(outpath, filename, "_geneInfo.csv"))
    write.table(output, file = paste0(filename, "_WGCNA.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
    save(data, softPower, MEs, dissTOM, moduleLabels, moduleColors, dynamicColors, geneTree, file = paste(outpath, filename, "_NetworkConstruction_ds-cutheight", ds,"-", cutHeight, ".RData", sep = ""))
}

# enrichR GO for modules
library(enrichR)
library(WriteXLS)

dbs <- c("KEGG_2021_Human", "GO_Biological_Process_2021", "BioCarta_2016", "Reactome_2022")

input_paths = "~/Desktop/HRM-THC-snRNA-seq/WGCNA-RNA/HGNC_outputs/"
files = list.files(input_paths)
outpath = "~/Desktop/HRM-THC-snRNA-seq/WGCNA-RNA/Module\ annotation/"

for (f in files){

    output <- read.delim(paste0(input_paths, f), stringsAsFactors = FALSE)
    filename <- gsub('_.*$', '', f)
    mods <- unique(output$MODULE)
    result = list()
    return.score <- ''

    count = 1
    for (n in mods){
        genes <- output$GENE[which(output$MODULE == n)]
        go.genes <- enrichr(genes, dbs)
        list.go.genes = do.call("rbind", go.genes)
        list.go.genes$padj <- p.adjust(list.go.genes$P.value, method = "BH")
        list.go.genes <- list.go.genes[order(list.go.genes$padj),]
        result[[n]] <- list.go.genes
        if (length(which(list.go.genes$padj < 0.05)) != 0){
            return.score[count] = 1
        }else{return.score[count] = 0}
        count = count + 1
    }
    return.score <- as.numeric(return.score)
    score <- round(sum(return.score)/length(return.score), 3)
    WriteXLS(result, paste0(outpath, filename, "_", score, "_enrichR_results.xlsx"), names(result), row.names = TRUE)
}
