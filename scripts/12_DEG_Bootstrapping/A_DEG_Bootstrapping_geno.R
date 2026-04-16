# run in parallel using seq 20 | xargs -n 1 -P 15 Rscript bootstrap_genotype.r


library(Seurat)
library(patchwork)
library(ggplot2)


load("20250324.RData")

args <- commandArgs(trailingOnly = TRUE)

argument1 <- args[1]
run_id <- args[1]
int_run_id <- as.integer(run_id)



objs <- ls()
seurat_obj_name <- objs[sapply(objs, function(x) inherits(get(x), "Seurat"))]
data <- get(seurat_obj_name)



DefaultAssay(data) <- "SCT"

counts <- GetAssayData(object = data, layer = "counts", assay = "RNA")

all_genes <- rownames(counts)

# 2. Create logical vectors for each gene type you want to remove
remove_rps <- startsWith(all_genes, "Rps")
remove_rpl <- startsWith(all_genes, "Rpl")
remove_mt  <- startsWith(all_genes, "mt-")

# 3. Combine them into a single vector.
# This vector is TRUE for any gene that matches ANY of the three patterns.
genes_to_remove <- remove_rps | remove_rpl | remove_mt

# 4. Now, subset your counts object ONE time using the inverted logical vector
# This keeps all rows where 'genes_to_remove' is FALSE.
counts_filtered <- counts[!genes_to_remove, ]
counts <- counts_filtered

data <- subset(data, features = rownames(counts))



sample_cells_by_type <- function(seurat_obj, column_name, value, n_cells) {
  matching_cells <- colnames(seurat_obj)[seurat_obj@meta.data[[column_name]] == value]
  sampled_cells <- sample(matching_cells, size = min(n_cells, length(matching_cells)),replace=TRUE)
  return(sampled_cells)
}

# this probably works, a bit manual
# duplicated cells names are dedupped using underscores, i.e. cellA, cellA, cellB ... -> cellA, cellA_1, cellB ...
generate_sample <- function(seurat_obj){
    inds <- append(sample_cells_by_type(seurat_obj, "Genotype", "Het", n_cells = 37747),
               sample_cells_by_type(seurat_obj, "Genotype", "WT", n_cells = 24164))
    counts.smp <- seurat_obj[["SCT"]]$counts[,inds]
    # even though only counts is used, data is required otherwise findmarkers errors
    data.smp <- seurat_obj[["SCT"]]$data[,inds]
    colnames(data.smp) <- make.unique(colnames(data.smp), sep="_")
    colnames(counts.smp) <-make.unique(colnames(counts.smp), sep="_")
    metadt <- data.frame(row.names = make.unique(inds, sep="_"),detailed.celltype = as.character(data$detailed.celltype[inds]),Genotype= data$Genotype[inds])
    return(CreateSeuratObject(data = data.smp, meta.data = metadt, counts=counts.smp))
}


OUTPUT_DIR <- "bootstrp_markers_geno/"

for (i in 1:50){
    # set the seed, very important
    # 100 * run id  + iteration: each iteration has unique seed
    set.seed(int_run_id*100 + i)
    print(paste0("Running ", i))
    sampled.obj <- generate_sample(data)
    Idents(object = sampled.obj) <- sampled.obj@meta.data[["detailed.celltype"]]
    genotype.markers <- list()
    for (clu in levels(Idents(sampled.obj))){
        # need try catch to deal w/ "biased" clusters, i.e. clusters w/ few Het (Ast-3)
        result <- tryCatch({
            # very important, min.pct, logfc threshold = 0 -> excluding no genes
            FindMarkers(object = sampled.obj, ident.1 = "Het", group.by = 'Genotype', subset.ident = clu,slot="counts",min.pct = 0, logfc.threshold = 0)
        }, error = function(e){
            NA
        })
        if (all(!is.na(result))){
            genotype.markers[[clu]] = result
        }

    }
    file_name <- paste0(OUTPUT_DIR, "geno_markers_run_id_",run_id, "_s_", i, ".rds")
    saveRDS(genotype.markers, file = file_name)
    rm(genotype.markers, sampled.obj)
    print(paste0("Finished ", i))

}
