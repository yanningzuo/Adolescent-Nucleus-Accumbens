set.seed(12)

# parallelization
library(future)
plan("multicore", workers = 8)

library(WGCNA)
library("flashClust")
options(stringsAsFactors = FALSE);
allowWGCNAThreads()

# reload data

input_paths = c("./WGCNA-RNA/Inputs/")
files = list.files(input_paths)

for (f in files){

    filename <- gsub('_.*$', '', f)
    setwd(paste0("./WGCNA-RNA/Data/", filename))
    files <- list.files(pattern = ".RData")
    load(files)

    # Quantifying module–trait associations
    # Define numbers of genes and samples
    nGenes = ncol(data)
    nSamples = nrow(data)

    # Recalculate MEs with color labels
    MEs0 = moduleEigengenes(data, moduleColors)$eigengenes
    MEs = orderMEs(MEs0)

    metadata <- read.csv("./WGCNA-SCT/metadata.csv")
    metadata <- metadata[which(metadata$X %in% rownames(data)), ]

    datTraits <- metadata[, c(3, 4)]
    moduleTraitCor = cor(MEs, datTraits, use = "p")
    moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
    write.csv(moduleTraitCor, file = paste0(filename, "_moduleTraitCor.csv"))
    write.csv(moduleTraitPvalue, file = paste0(filename, "_moduleTraitPvalue.csv"))

}

# ============================================================
# Visualization: module-trait correlation heatmaps + bar plot
# ============================================================

cor.sig <- read.csv("./WGCNA-RNA/all-sig_moduleTraitCor.csv")
mod.sig <- read.csv("./WGCNA-RNA/all-sig_moduleTraitPvalue.csv")
cor.sig <- cor.sig[, -1]; mod.sig <- mod.sig[, -1]

### heatmaps

library(pheatmap)
library(RColorBrewer)
library(reshape)

# https://stackoverflow.com/questions/15505607/diagonal-labels-orientation-on-x-axis-in-heatmaps

library(grid)     ## Need to attach (and not just load) grid package
## For pheatmap_1.0.8 and later:
draw_colnames_45 <- function (coln, gaps, ...) {
    coord = pheatmap:::find_coordinates(length(coln), gaps)
    x = coord$coord - 0.5 * coord$size
    res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 1.5, hjust = 1, rot = 45, gp = gpar(...))
    return(res)}

## 'Overwrite' default draw_colnames with your own version
## (value must be the function object, not a quoted string -- assignInNamespace()
## with a string here would make every later pheatmap() call crash, since
## pheatmap calls draw_colnames() internally to render column labels)
assignInNamespace(x="draw_colnames", value=draw_colnames_45,
ns=asNamespace("pheatmap"))

# Z score across all module's correlation r

input_paths = c("./WGCNA-RNA/Data/")
files = list.files(input_paths)

mod <- data.frame()
cor <- data.frame()

for (f in files){

    pv.file <- read.csv(paste0(input_paths, f, "/", f, "_moduleTraitPvalue.csv"))
    cor.file <- read.csv(paste0(input_paths, f, "/", f, "_moduleTraitCor.csv"))

    if (length(which(pv.file$X == "MEgrey")) != 0){
        pv.file <- pv.file[-which(pv.file$X == "MEgrey"), ]
        cor.file <- cor.file[-which(cor.file$X == "MEgrey"), ]
    }

    pv.file$X <- gsub("ME", "", pv.file$X)
    pv.file$X <- paste0(f, "_", pv.file$X)

    cor.file$X <- gsub("ME", "", cor.file$X)
    cor.file$X <- paste0(f, "_", cor.file$X)

    mod <- rbind(mod, pv.file)
    cor <- rbind(cor, cor.file)

}

write.csv(mod, "./WGCNA-RNA/all-moduleTraitPvalue.csv")
write.csv(cor, "./WGCNA-RNA/all-moduleTraitCor.csv")

z.all = data.frame(Genotype = (cor$Genotype - mean(cor$Genotype))/sd(cor$Genotype), Treatment = (cor$Treatment - mean(cor$Treatment))/sd(cor$Treatment))
rownames(z.all) = cor$X
write.csv(mod, "./WGCNA-RNA/all-moduleCor-Z.csv")

# all sig modules

sig.mod = mod$X[which(mod$Genotype < 0.05 | mod$Treatment < 0.05)] # 43
sig.cor = data.frame(Mod = sig.mod, Treat.Cor = cor$Treatment[which(cor$X %in% sig.mod)], Geno.Cor = cor$Genotype[which(cor$X %in% sig.mod)])
sig.Z = data.frame(Treat.Z = z.all$Treatment[which(rownames(z.all) %in% sig.mod)], Geno.Z = z.all$Genotype[which(rownames(z.all) %in% sig.mod)])
rownames(sig.Z) = sig.mod
sig.p = data.frame(Mod = sig.mod, Treat.P = mod$Treatment[which(mod$X %in% sig.mod)], Geno.P = mod$Genotype[which(mod$X %in% sig.mod)])

test_labels <- sig.p
rownames(test_labels) <- sig.p[, 1]
test_labels <- test_labels[, -1]
test_labels <- t(test_labels)

test_labels <- as.matrix(test_labels)
test_labels[test_labels < 0.001] <- "***"
test_labels[test_labels < 0.01 & test_labels >= 0.001] <- "**"
test_labels[test_labels < 0.05 & test_labels >= 0.01] <- "*"
test_labels[test_labels >= 0.05] <- ""

breaksList = seq(-3, 3, by = 1)
pheatmap(t(sig.Z), scale = "none", trace = "none", col = colorRampPalette(rev(brewer.pal(9, "RdBu")))(length(breaksList)), breaks = breaksList, border_color = NA, display_numbers = test_labels, fontsize_number = 12, cellheight = 15, cellwidth = 15)

# treatment + genotype

sig.mod = mod$X[which(mod$Genotype < 0.05 & mod$Treatment < 0.05)]
sig.cor = data.frame(Mod = sig.mod, Treat.Cor = cor$Treatment[which(cor$X %in% sig.mod)], Geno.Cor = cor$Genotype[which(cor$X %in% sig.mod)])
sig.Z = data.frame(Treat.Z = z.all$Treatment[which(rownames(z.all) %in% sig.mod)], Geno.Z = z.all$Genotype[which(rownames(z.all) %in% sig.mod)])
rownames(sig.Z) = sig.mod
sig.p = data.frame(Mod = sig.mod, Treat.P = mod$Treatment[which(mod$X %in% sig.mod)], Geno.P = mod$Genotype[which(mod$X %in% sig.mod)])

test_labels <- sig.p
rownames(test_labels) <- sig.p[, 1]
test_labels <- test_labels[, -1]
test_labels <- t(test_labels)

test_labels <- as.matrix(test_labels)
test_labels[test_labels < 0.001] <- "***"
test_labels[test_labels < 0.01 & test_labels >= 0.001] <- "**"
test_labels[test_labels < 0.05 & test_labels >= 0.01] <- "*"
test_labels[test_labels >= 0.05] <- ""

breaksList = seq(-3, 3, by = 1)
pheatmap(t(sig.Z), scale = "none", trace = "none", col = colorRampPalette(rev(brewer.pal(9, "RdBu")))(length(breaksList)), breaks = breaksList, border_color = NA, display_numbers = test_labels, fontsize_number = 16, cellheight = 32, cellwidth = 32)

## bar plot of trait correlation

library(ggplot2)
library(stringr)

mod.sig$Cell.type = str_split_i(mod.sig$X, "_", 1)
Cell.type <- rep(unique(mod.sig$Cell.type), 2)
Correlation <- c(rep("Treatment", length(unique(mod.sig$Cell.type))), rep("Genotype", length(unique(mod.sig$Cell.type))))
Count <- rep(0, length(unique(mod.sig$Cell.type))*2)

df <- data.frame(Cell.type, Correlation, Count)
df$Cell.type <- factor(df$Cell.type, levels = c("D1.MSN-1", "D1.MSN-2", "D1.MSN-4", "D1.MSN-5", "D2.MSN-1", "D3.MSN", "GABA-1", "GABA-2", "GABA-3", "GABA-4", "Ast-3", "Ol-2"))

for (i in unique(mod.sig$Cell.type)){
    items = which(mod.sig$Cell.type == i)
    temp = mod.sig[items, ]
    df$Count[which(df$Correlation == "Treatment" & df$Cell.type == i)] = length(which(temp$Treatment < 0.05))
    df$Count[which(df$Correlation == "Genotype" & df$Cell.type == i)] = length(which(temp$Genotype < 0.05))
}

ggplot(df, aes(fill = Correlation, y = Count, x = Cell.type)) +
    geom_bar(position = "dodge", stat = "identity") +
    theme_minimal() +
    labs(x= "", y = "Module counts") +
    scale_fill_manual(values = c("#00bfc4", "#f8766d")) +
    theme(text = element_text(size = 16), axis.text.x = element_text(size = 12, angle = 45, hjust=1)) +
    ggtitle("Trait-module correlation")

# ============================================================
# Visualization: module eigengene (ME) clustering tree figures
# ============================================================

library(stringr)
library(ggplotify)
library(ggpubr)

input.path = "./WGCNA-RNA/Data/"
files = list.files(input.path, recursive = TRUE, pattern = ".RData")
files = files[-c(20, 23)]

pics <- list()

for (f in files){

    load(paste0(input.path, f))
    filename = str_split_i(f, "/", 1)

    MEDiss = 1 - cor(MEs)
    METree = hclust(as.dist(MEDiss), method = "average")

    plot(METree, main = filename, xlab = "", sub = "")
    pics[[filename]] = recordPlot()

    rm(softPower, dynamicMods, MEList, MEs, MEDiss, METree, dissTOM, merge, moduleLabels, moduleColors, dynamicColors, geneTree, geneModulebicor, geneModuleMembership, MMPvalue, geneInfo0, geneInfo, output, data)
    gc(TRUE)

}

ggpubr::ggarrange(plotlist = pics, nrow = 9, ncol = 3, byrow = TRUE) # not working, possibly big size

pics <- list()

for (f in files){

    load(paste0(input.path, f))
    filename = str_split_i(f, "/", 1)

    MEDiss = 1 - cor(MEs)
    METree = hclust(as.dist(MEDiss), method = "average")

    pdf(file = paste0("./WGCNA-RNA/MEtree/", filename, ".pdf"), width = 20.31, height = 5.73)
    plot(METree, main = filename, xlab = "", sub = "")
    dev.off()

    rm(softPower, dynamicMods, MEList, MEs, MEDiss, METree, dissTOM, merge, moduleLabels, moduleColors, dynamicColors, geneTree, geneModulebicor, geneModuleMembership, MMPvalue, geneInfo0, geneInfo, output, data)
    gc(TRUE)

}

