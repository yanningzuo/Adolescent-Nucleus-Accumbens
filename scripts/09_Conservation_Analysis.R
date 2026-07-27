###############################################################################
# 09 — CONSERVATION ANALYSES
#
#   PART 1 — CROSS-AGE:     adolescent mouse NAc (this study) vs.
#                           adult mouse NAc (GSE118020, P70)
#
#   PART 2 — CROSS-SPECIES: adolescent mouse NAc (this study) vs.
#                           adult human NAc (Tran et al., n = 8 donors)
#                           -> manuscript Fig. S10
#
# Manuscript: Zuo et al., bioRxiv 2025.04.03.646846
# Methods, "Cross-species conservation analysis":
#   "We calculated cell type-specific t-statistics on SCTransformed expression
#    values following established methods... We used the cell type identities
#    reported in the original study and applied SCTransform normalization,
#    regressing out sequencing depth, donor, sex, processing batch, and
#    ribosomal gene content (percent.RPL and percent.RPS). We calculated
#    Pearson correlations using the top 100 upregulated genes per cell type."
#
# Both parts use the IDENTICAL statistical approach — same t-statistic, same
# top-100-marker selection, same pairwise Pearson correlation, same heatmap.
# They share the helper functions defined once below. Part 2 adds exactly two
# things Part 1 does not need, because Part 1 is mouse-vs-mouse:
#   (i)  the human dataset is re-processed from raw counts (QC + SCTransform
#        with the donor / sex / batch covariates a multi-donor human dataset
#        requires), and
#   (ii) an ortholog-mapping step converts human gene symbols to mouse symbols
#        before the two matrices can be intersected.
#
# Part 1 and Part 2 are independent and can be run separately, EXCEPT that
# Part 2 reuses the adolescent-mouse t-statistic matrix produced in Part 1
# step 1a. Run that step first if you only want Part 2.
#
#------------------------------------------------------------------------------
# METHOD OVERVIEW (three steps, both parts)
#
#   STEP 1. For each dataset independently, compute a cell type-specificity
#           t-statistic for every gene x cell type (one-vs-rest).
#   STEP 2. Restrict both datasets to a shared gene set: genes present in both,
#           further restricted to the union of the top 100 marker genes per
#           adolescent NAc cell type. (Part 2 first maps human symbols to mouse.)
#   STEP 3. Pearson-correlate every mouse cell type's t-statistic vector against
#           every other-dataset cell type's vector; plot as a heatmap.
#
#------------------------------------------------------------------------------
# MARKER GENE SELECTION — answers to the questions raised in review
#
#   Q: How was the marker gene list selected?
#   A: The marker genes are NOT computed in this script. They are read from the
#      `markers` list object created in 03_SCT_Clustering_Annotation.R:
#
#          DefaultAssay(data) <- "SCT"          # <- set before the loop
#          markers <- list()
#          for (clu in levels(Idents(data))){
#              markers[[clu]] <- FindMarkers(object = data, ident.1 = clu)
#          }
#
#      i.e. Seurat FindMarkers(), one cluster vs. all remaining cells, run once
#      per cell type and stored in a named list. That list is saved inside the
#      annotated .RData snapshot and is simply re-loaded here — the conservation
#      analysis does NOT re-run FindMarkers().
#
#   Q: Which assay was FindMarkers() run on? Was it SCTransformed data?
#   A: The SCT assay. DefaultAssay(data) was set to "SCT" before the marker loop
#      (03_SCT_Clustering_Annotation.R, line 19; loop at lines 46-49), and
#      FindMarkers() was called without an `assay =` argument, so it used the
#      default assay = SCT. All other FindMarkers() parameters were left at
#      Seurat defaults:
#        test.use        = "wilcox"   (Wilcoxon rank sum)
#        logfc.threshold = 0.25
#        min.pct         = 0.1
#      Note this is the *marker* (one-vs-rest) call. It is a different call from
#      the treatment/genotype DEG analysis in 04_DEG_Pathway_Analysis.R.
#
#   Q: Were the top 100 genes selected based on log fold change?
#   A: Yes. Within each cell type the marker table is (1) filtered to
#      p_val_adj < 0.05 AND avg_log2FC > 0 (upregulated markers only), then
#      (2) sorted by descending avg_log2FC, then (3) the top 100 taken.
#      Selection is by avg_log2FC, not by p-value and not by pct.diff.
#
#      The per-cell-type lists are concatenated into ONE vector and
#      deduplicated, so the final gene set is the *union* of all adolescent NAc cell types'
#      top-100 lists — not 100 genes total, and not a separate gene set per
#      correlation. Every cell of the correlation matrix uses this same set.
#
#      The significance filter is applied BEFORE the top-100 cut, so the 100 are
#      drawn only from genes already passing p_val_adj < 0.05 & avg_log2FC > 0.
#      A cell type with fewer than 100 such genes contributes all of them and
#      nothing more — no additional genes are pulled from further down the
#      unfiltered marker table to top the list back up to 100.
#
#      In PART 2 the markers are taken from the MOUSE side only. Human marker
#      genes are never computed; the human contribution is its t-statistic
#      matrix, ortholog-mapped and then subset to the mouse-derived gene set.
#
###############################################################################

library(Seurat)
library(SeuratObject)
library(dplyr)
library(ggplot2)
library(stringr)
library(patchwork)
library(pheatmap)
library(RColorBrewer)

set.seed(12)

# parallelization
library(future)
plan("multicore", workers = 8)
plan()
options(future.globals.maxSize = 480 * 1024 ^ 3)   # for 480 Gb RAM


###############################################################################
# INPUT / OUTPUT PATHS
#
# All paths are relative to the project root, so run this script with the
# project root as the working directory. (Deliberately no setwd() calls — the
# earlier versions of this script mixed setwd() with root-relative paths, which
# silently broke the second and third load() calls.)
###############################################################################

# --- adolescent mouse NAc (this study) ---------------------------------------
# Annotated Seurat object from 03_SCT_Clustering_Annotation.R. Supplies both the
# expression data for the t-statistics and the `markers` list used in step 2.
mouse_obj_crossage    <- "./Merged/20230123_SCT_noGlut-61911.RData"

# Snapshot used for the cross-species marker lookup. This is a LATER snapshot
# whose cell type labels differ slightly (see the harmonization note in Part 2).
mouse_obj_crossspecies <- "./Merged/20250130.RData"

mouse_tstats_file     <- "./Merged/t-stats.Rdata"

# --- adult mouse NAc (GSE118020, P70) ----------------------------------------
adult_expr_file       <- "./Adult-NAc/GSE118020_expression_matrix.csv.gz"
adult_meta_file       <- "./Adult-NAc/GSE118020_meta_data.csv.gz"
adult_tstats_file     <- "./Adult-NAc/adult_t-stats.Rdata"

# --- adult human NAc (Tran et al.) -------------------------------------------
# Published as a SingleCellExperiment. Place the .rda here before running Part 2.
human_sce_file        <- "./Human-NAc/SCE_NAc-n8_tran-etal.rda"
human_tstats_file     <- "./Human-NAc/t-stats.Rdata"


###############################################################################
# SHARED HELPER FUNCTIONS
###############################################################################

#' Cell type-specificity t-statistic
#'
#' For each gene and each cell type, compares expression in that cell type
#' against all other cells (one-vs-rest):
#'
#'     s_pooled = sqrt( (sd_celltype^2 + sd_other^2) / 2 )
#'     t        = sqrt(n_celltype) * (mean_celltype - mean_other) / s_pooled
#'
#' High positive t = gene is specifically enriched in that cell type. The
#' resulting matrix (genes x cell types) is the input to the correlation step.
#'
#' @param seurat_obj Seurat object with Idents() already set to the cell type
#'                   labels you want as columns.
#' @param assay      Assay whose @data (normalized) slot is used. "SCT"
#'                   throughout this study.
#' @return data.frame, rows = genes, cols = cell types.
celltype_tstats <- function(seurat_obj, assay = "SCT") {

    metadata <- FetchData(seurat_obj, vars = "ident")
    exp      <- seurat_obj@assays[[assay]]@data

    t.stats <- data.frame(Genes = rownames(exp))

    for (ident in unique(seurat_obj@active.ident)) {

        cells    <- rownames(metadata)[which(metadata$ident == ident)]
        temp.exp <- exp[, which(colnames(exp) %in% cells)]

        other.cells <- rownames(metadata)[which(metadata$ident != ident)]
        other.exp   <- exp[, which(colnames(exp) %in% other.cells)]

        s.celltype <- apply(temp.exp,  1, sd)
        s.other    <- apply(other.exp, 1, sd)
        s.pooled   <- sqrt((s.celltype^2 + s.other^2) / 2)

        m.celltype <- apply(temp.exp,  1, mean)
        m.other    <- apply(other.exp, 1, mean)

        t <- sqrt(length(cells)) * (m.celltype - m.other) / s.pooled

        t <- as.data.frame(t)
        colnames(t) <- ident
        t.stats <- cbind(t.stats, t)
    }

    t.stats <- t.stats[, -1]   # drop the Genes seed column
    return(t.stats)
}


#' Union of the top-N upregulated marker genes across all cell types
#'
#' See the "MARKER GENE SELECTION" note in the header for full provenance.
#'
#' @param markers_list Named list of FindMarkers() output, one element per cell
#'                     type (created in 03_SCT_Clustering_Annotation.R, on the
#'                     SCT assay, one-vs-rest, Seurat defaults).
#' @param cell_types   Which cell types to draw markers from. These must all be
#'                     present in names(markers_list).
#' @param n            Number of top genes per cell type. 100 in this study.
#' @return character vector of unique gene symbols.
top_markers_union <- function(markers_list, cell_types, n = 100) {

    missing <- setdiff(cell_types, names(markers_list))
    if (length(missing) > 0) {
        stop("No marker table for cell type(s): ", paste(missing, collapse = ", "),
             ". The markers snapshot and the t-statistic matrix must use the ",
             "same cell type labels.")
    }

    marker <- character(0)

    for (ident in cell_types) {
        temp <- markers_list[[ident]]

        # upregulated markers only, significant after FDR correction
        temp <- temp[which(temp$p_val_adj < 0.05 & temp$avg_log2FC > 0), ]

        # rank by effect size, take the top n
        temp  <- temp[order(-temp$avg_log2FC), ]
        # head() rather than [1:n]: cell types with fewer than n genes surviving
        # the filter above contribute all of those and nothing more, instead of
        # being NA-padded out to length n
        genes <- head(rownames(temp), n)

        marker <- c(marker, genes)
    }

    return(unique(marker))
}


#' Restrict two t-stat matrices to shared genes and put them in identical order
#'
#' This step is essential. cor.test() pairs values by POSITION, not by name, so
#' subsetting to shared genes alone is NOT sufficient — the two datasets do not
#' store genes in the same order, and correlating them unsorted silently
#' compares mismatched gene pairs.
align_by_gene <- function(a, b) {

    shared.genes <- intersect(rownames(a), rownames(b))

    a <- a[which(rownames(a) %in% shared.genes), ]
    b <- b[which(rownames(b) %in% shared.genes), ]

    a <- a[order(rownames(a)), ]
    b <- b[order(rownames(b)), ]

    stopifnot(identical(rownames(a), rownames(b)))

    return(list(a = a, b = b))
}


#' Pairwise Pearson correlation between two t-statistic matrices
#'
#' Inputs must already be gene-aligned (see align_by_gene()).
#'
#' @return list(cor = matrix, p = matrix); rows = cell types of `b`,
#'         cols = cell types of `a`.
pairwise_cor <- function(a, b) {

    cor.mat <- matrix(0, nrow = ncol(b), ncol = ncol(a),
                      dimnames = list(colnames(b), colnames(a)))
    p.mat   <- cor.mat

    for (i in 1:ncol(a)) {
        for (j in 1:ncol(b)) {
            res <- cor.test(a[, i], b[, j], method = "pearson")
            cor.mat[j, i] <- res[["estimate"]]
            p.mat[j, i]   <- res[["p.value"]]
        }
    }

    return(list(cor = cor.mat, p = p.mat))
}


#' Conservation heatmap
#'
#' Cells are labeled with the correlation coefficient; non-significant cells
#' (p >= 0.05) are left blank.
conservation_heatmap <- function(cor.mat, p.mat, fontsize_number = 12,
                                 cellheight = NA, cellwidth = NA) {

    test_labels <- round(cor.mat, 2)
    test_labels[which(p.mat >= 0.05)] <- " "

    breaksList <- seq(-0.5, 0.5, by = 0.01)

    pheatmap(cor.mat,
             scale           = "none",
             trace           = "none",
             breaks          = breaksList,
             col             = colorRampPalette(rev(brewer.pal(9, "RdBu")))(length(breaksList)),
             cluster_cols    = TRUE,
             cluster_rows    = TRUE,
             border_color    = NA,
             display_numbers = test_labels,
             fontsize_number = fontsize_number,
             cellheight      = cellheight,
             cellwidth       = cellwidth)
}


#' Rotate pheatmap column labels 45 degrees (used for the wider Part 2 heatmap)
#' https://stackoverflow.com/questions/15505607/diagonal-labels-orientation-on-x-axis-in-heatmaps
rotate_pheatmap_colnames <- function() {
    library(grid)   ## Need to attach (and not just load) grid package
    draw_colnames_45 <- function(coln, gaps, ...) {
        coord <- pheatmap:::find_coordinates(length(coln), gaps)
        x <- coord$coord - 0.5 * coord$size
        res <- textGrob(coln, x = x, y = unit(1, "npc") - unit(3, "bigpts"),
                        vjust = 1.5, hjust = 1, rot = 45, gp = gpar(...))
        return(res)
    }
    assignInNamespace(x = "draw_colnames", value = draw_colnames_45,
                      ns = asNamespace("pheatmap"))
}


#' Map human gene symbols to mouse gene symbols via MGI homology classes
#'
#' Uses the MGI HOM_MouseHumanSequence report: two symbols are orthologs if they
#' share a DB.Class.Key. Symbols with no mouse ortholog are returned unchanged,
#' so they simply fail to intersect with the mouse matrix later and drop out.
#'
#' NOTE: this fetches a live file from informatics.jax.org over plain http. If
#' that URL or file format changes upstream, results could shift. To make a run
#' fully reproducible, download the report once and read it from disk instead.
convert_human_to_mouse <- function(gene_list) {

    require("dplyr")
    mouse_human_genes <- read.csv(
        "http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",
        sep = "\t"
    )

    output <- data.frame(Input = gene_list, Output = rep("", length(gene_list)))

    for (gene in gene_list) {
        class_key <- (mouse_human_genes %>%
                      filter(Symbol == gene & Common.Organism.Name == "human"))[['DB.Class.Key']]
        if (!identical(class_key, integer(0))) {
            mouse_genes <- (mouse_human_genes %>%
                            filter(DB.Class.Key == class_key &
                                   Common.Organism.Name == "mouse, laboratory"))[, "Symbol"]
            for (mouse_gene in mouse_genes) {
                output$Output[which(output$Input == gene)] <- mouse_gene
            }
        }
    }

    no.gene <- which(output$Output == "")
    output$Output[no.gene] <- output$Input[no.gene]
    return(output)
}



###############################################################################
###############################################################################
##                                                                           ##
##   PART 1 — CROSS-AGE CONSERVATION                                         ##
##            adolescent mouse NAc  vs.  adult mouse NAc (GSE118020)         ##
##                                                                           ##
###############################################################################
###############################################################################


#==============================================================================
# 1.1 — t-statistics for the adolescent mouse NAc (this study)
#
# Also reused by PART 2, so run this step even if you only want cross-species.
#==============================================================================

load(mouse_obj_crossage)   # provides `data` (Seurat object) and `markers` (list)

DefaultAssay(data)    <- "SCT"
Idents(object = data) <- data@meta.data$detailed.celltype

t.stats <- celltype_tstats(data, assay = "SCT")
save(t.stats, file = mouse_tstats_file)


#==============================================================================
# 1.2 — t-statistics for the adult mouse NAc (GSE118020, P70)
#==============================================================================

data <- read.csv(adult_expr_file)
rownames(data) <- data$X
data <- data[, -1]
metadata <- read.csv(adult_meta_file)

data <- CreateSeuratObject(counts = data, project = "Adult-NAc")
data <- AddMetaData(
    object   = data,
    metadata = data.frame(detailed.ident = metadata$L2_Clustering,
                          row.names = rownames(data@meta.data))
)

# Normalized the same way as our own data so the two t-stat matrices are
# comparable. Only sequencing depth is regressed out here, as this public
# dataset does not carry the donor/batch covariates our object has.
data <- SCTransform(data, vars.to.regress = c("nCount_RNA"), verbose = TRUE)

DefaultAssay(data)    <- "SCT"
Idents(object = data) <- data@meta.data[["detailed.ident"]]

t.stats <- celltype_tstats(data, assay = "SCT")
save(t.stats, file = adult_tstats_file)


#==============================================================================
# 1.3 — build the shared gene set
#==============================================================================

load(mouse_tstats_file)
adolescence <- t.stats

load(adult_tstats_file)
adult <- t.stats
rm(t.stats)

# (a) genes present in both datasets, in identical row order
aligned     <- align_by_gene(adolescence, adult)
adolescence <- aligned$a
adult       <- aligned$b

# (b) restrict to the union of the top 100 upregulated markers per adolescent
#     mouse cell type. `markers` comes from the snapshot loaded in step 1.1 and
#     was produced by FindMarkers() on the SCT assay — see header for details.
load(mouse_obj_crossage)   # reload: `markers` was overwritten above

marker <- top_markers_union(markers,
                            cell_types = colnames(adolescence),
                            n = 100)

adolescence <- adolescence[which(rownames(adolescence) %in% marker), ]
adult       <- adult[which(rownames(adult) %in% marker), ]


#==============================================================================
# 1.4 — correlate every adolescent cell type against every adult cell type
#==============================================================================

res <- pairwise_cor(adolescence, adult)

conservation_heatmap(res$cor, res$p, fontsize_number = 8)



###############################################################################
###############################################################################
##                                                                           ##
##   PART 2 — CROSS-SPECIES CONSERVATION            (manuscript Fig. S10)    ##
##            adolescent mouse NAc  vs.  adult human NAc (Tran et al.)       ##
##                                                                           ##
###############################################################################
###############################################################################


#==============================================================================
# 2.1 — t-statistics for the adult human NAc (Tran et al.)
#==============================================================================

# The published dataset ships as a SingleCellExperiment; convert to Seurat.
# After conversion the counts live in an assay named "originalexp", and the
# authors' cell type labels are in meta.data$cellType.
load(human_sce_file)
data <- as.Seurat(sce.nac.tran)

Idents(object = data) <- data@meta.data[["cellType"]]

# --- QC ---------------------------------------------------------------------
# Filters mirror those applied to our own mouse data in 02_QC_Doublet_Removal.R,
# with thresholds adjusted for this dataset's depth. Note human gene symbols are
# uppercase, hence "^MT-" / "^RPL" / "^RPS" rather than the mouse patterns.
data[["percent.MT"]] <- PercentageFeatureSet(data, pattern = "^MT-")
data <- PercentageFeatureSet(data, pattern = "^RPL", col.name = "percent.RPL")
data <- PercentageFeatureSet(data, pattern = "^RPS", col.name = "percent.RPS")

data <- subset(data, subset = nFeature_originalexp > 500
                            & nFeature_originalexp < 8000
                            & nCount_originalexp   > 1000
                            & nCount_originalexp   < 30000
                            & percent.RPS < 1
                            & percent.RPL < 1
                            & percent.MT  < 1)

# --- Normalization ----------------------------------------------------------
# Multi-donor human dataset, so donor / sex / processing batch are regressed out
# in addition to sequencing depth and ribosomal content. (The Methods text lists
# depth, donor, sex, batch and ribosomal content; percent.MT is also regressed
# here, as in the original run.)
data <- SCTransform(
    data,
    assay           = "originalexp",
    vars.to.regress = c("donor", "sex", "processBatch",
                        "percent.MT", "nCount_originalexp",
                        "percent.RPL", "percent.RPS"),
    verbose = TRUE
)

data <- RunPCA(data, verbose = TRUE)
data <- FindNeighbors(data, dims = 1:50, verbose = TRUE)
data <- FindClusters(data, verbose = TRUE, resolution = 0.7)
data <- RunUMAP(data, dims = 1:50, verbose = TRUE)

# Use the ORIGINAL STUDY's cell type labels for the comparison (per Methods),
# not our re-clustering. The re-clustering above is only run so the object
# carries a PCA/UMAP; the correlation uses `cellType`.
Idents(object = data) <- data@meta.data[["cellType"]]

# Drop the classes the original authors flagged as low quality / doublets.
data <- subset(data,
               idents = c("drop.doublet_A", "drop.doublet_B", "drop.doublet_C",
                          "drop.doublet_D", "drop.lowNTx"),
               invert = TRUE)

t.stats <- celltype_tstats(data, assay = "SCT")
save(t.stats, file = human_tstats_file)

human <- t.stats


#==============================================================================
# 2.2 — load the adolescent mouse t-statistics and marker list
#==============================================================================

# Mouse t-statistics computed in step 1.1 — same matrix, no need to recompute.
load(mouse_tstats_file)
mouse <- t.stats
rm(t.stats)

# Marker list for the cross-species comparison. NOTE this is a LATER annotation
# snapshot than the one used in Part 1, so its cell type labels differ slightly
# from the column names of the mouse t-statistic matrix — harmonized in 2.3(b).
load(mouse_obj_crossspecies)   # provides `markers`


#==============================================================================
# 2.3 — map human genes to mouse orthologs, then build the shared gene set
#==============================================================================

human.genes <- convert_human_to_mouse(rownames(human))$Output

# Several human genes can map to the same mouse ortholog. Where that happens,
# the t-statistics of the duplicate rows are SUMMED and collapsed into a single
# row, so that each mouse symbol appears exactly once.
#
# Caveat worth knowing: summing t-statistics is not a variance-weighted
# combination, so a mouse gene backed by several human paralogs gets a larger
# magnitude than one backed by a single ortholog. Affected genes are a small
# minority of the matrix, but if you revisit this analysis, taking the mean (or
# the max-|t|) would be a defensible alternative.
dups      <- which(duplicated(human.genes))
dup.genes <- human.genes[duplicated(human.genes)]

for (gene in dup.genes) {
    nrow    <- which(human.genes == gene)
    add.num <- apply(human[nrow, ], 2, sum)
    human[nrow, ] <- add.num          # write the group sum into every member row
}

human           <- human[-dups, ]                   # keep the first of each group
human.genes     <- human.genes[!duplicated(human.genes)]
rownames(human) <- human.genes

# (a) genes present in both species, in identical row order
aligned <- align_by_gene(mouse, human)
mouse   <- aligned$a
human   <- aligned$b

# (b) harmonize the one cell type label that differs between the mouse t-stat
#     matrix (computed from the earlier snapshot, where this cluster is "NB")
#     and the markers list loaded in 2.2 (where it is "NPC"), so that the marker
#     lookup below resolves for every column
colnames(mouse)[which(colnames(mouse) == "NB")] <- "NPC"
names(markers)[which(names(markers)   == "NB")] <- "NPC"

# (c) restrict to the union of the top 100 upregulated markers per MOUSE cell
#     type — see header for how these were generated
marker <- top_markers_union(markers,
                            cell_types = colnames(mouse),
                            n = 100)

mouse <- mouse[which(rownames(mouse) %in% marker), ]
human <- human[which(rownames(human) %in% marker), ]


#==============================================================================
# 2.4 — correlate every mouse cell type against every human cell type
#==============================================================================

res <- pairwise_cor(mouse, human)

rotate_pheatmap_colnames()

conservation_heatmap(res$cor, res$p,
                     fontsize_number = 12,
                     cellheight = 30, cellwidth = 30)
