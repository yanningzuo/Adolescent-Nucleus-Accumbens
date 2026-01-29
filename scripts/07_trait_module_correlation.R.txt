
set.seed(12)

# parallelization
library(future)
plan("multicore", workers = 8)

library(WGCNA)
library("flashClust")
options(stringsAsFactors = FALSE);
allowWGCNAThreads() 

# reload data

input_paths = c("~/Desktop/HRM-THC-snRNA-seq/WGCNA-RNA/Inputs/")
files = list.files(input_paths)

for (f in files){

    filename <- gsub('_.*$', '', f)
    setwd(paste0("~/Desktop/HRM-THC-snRNA-seq/WGCNA-RNA/Data/", filename))
    files <- list.files(pattern = ".RData")
    load(files)

    # Quantifying module–trait associations
    # Define numbers of genes and samples
    nGenes = ncol(data)
    nSamples = nrow(data)

    # Recalculate MEs with color labels
    MEs0 = moduleEigengenes(data, moduleColors)$eigengenes
    MEs = orderMEs(MEs0)

    metadata <- read.csv("~/Desktop/HRM-THC-snRNA-seq/WGCNA-SCT/metadata.csv")
    metadata <- metadata[which(metadata$X %in% rownames(data)), ]

    datTraits <- metadata[, c(3, 4)]
    moduleTraitCor = cor(MEs, datTraits, use = "p")
    moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
    write.csv(moduleTraitCor, file = paste0(filename, "_moduleTraitCor.csv"))
    write.csv(moduleTraitPvalue, file = paste0(filename, "_moduleTraitPvalue.csv"))

}
