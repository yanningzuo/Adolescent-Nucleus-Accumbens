### Mergeomics using OUD as an example
## Marker set filtering (MDF)

# prep GWAS summary stats

GWAS <- read.delim("opi.DEPvEXP_trans.ancestry.noAF.tbl", stringsAsFactors = FALSE)
GWAS = GWAS[, c("rsID","P.value")]
GWAS$P.value = -log10(GWAS$P.value)
colnames(GWAS) = c("LOCUS", "VALUE")
write.table(GWAS, "2020_OUD.txt", row.names = FALSE, quote = FALSE, sep = "\t")

# MDF

cd ./Mergeomics/OUD/MDF/bash
chmod +x ./Mergeomics/bin/ldprune
export PATH=$PATH:./Mergeomics/bin/

export PATH=$PATH:./Mergeomics/OUD/MDF/bash/
chmod u+x oud_distance_20kb.bash
./Mergeomics/OUD/MDF/bash/oud_distance_20kb.bash

chmod u+x oud_eqtl-nac.bash
oud_eqtl-nac.bash

chmod u+x oud_sqtl-nac.bash
oud_sqtl-nac.bash

## Marker set enrichment analysis (MSEA)

library(ggtern)
source("./Mergeomics/Mergeomics_YZ_2020_v2.R")
path = "./Mergeomics/OUD/"
trait.item = 'OUD'
trim = '0.005'
trim = as.numeric(trim)
name = 'NAc-module-DEG'
dir.create(file.path(path, name))
setwd(file.path(path, name))

mapping = c("distance", "NAc-eqtl", "NAc-sqtl")
for(item in mapping){
    job.ssea <- list()
    job.ssea$label <- paste(trait.item, item, sep = ".")
    job.ssea$folder <- file.path(path, name)
    job.ssea$genfile <- paste0(path, "MDF/", item, "_genes.txt")
    job.ssea$locfile <- paste0(path, "MDF/", item, "_loci.txt")
    job.ssea$modfile <- paste0("./Mergeomics/", name, "_combined_module.txt")
    job.ssea$inffile <- paste0("./Mergeomics/", name, "_combined_module_info.txt")
    job.ssea$permtype <- "gene"	#optional
    job.ssea$maxgenes <- 1000	#optional
    job.ssea$nperm <- 10000	#optional
    job.ssea <- ssea.start(job.ssea)
    job.ssea <- ssea.prepare(job.ssea)
    job.ssea <- ssea.control(job.ssea)
    job.ssea <- ssea.analyze(job.ssea, trim_start = trim, trim_end =(1 - trim))
    job.ssea <- ssea.finish(job.ssea)
}

## post-MSEA analysis
# read in results files from MSEA
results = data.frame(stringsAsFactors = FALSE)
for (item in mapping){
    temp = read.delim(paste0(path, name, "/msea/", trait.item, ".", item, ".results.txt"), stringsAsFactors = FALSE)
    temp$Mapping = item
    results <- rbind(results, temp)
}

results = results[results$FDR<0.05 & results$MODULE!="_ctrlA" & results$MODULE!="_ctrlB",]
write.table(results, "MSEA_compiled_results.txt", row.names = FALSE, quote = FALSE, sep = "\t")

# trim reoccuring modules
reoccuring = results$MODULE[duplicated(results$MODULE)]
for(row in 1:nrow(results)){
  if(sum(results$MODULE[row]==reoccuring)==0){
    next
  }
  else{
    mappings = results$Mapping[results$MODULE==results$MODULE[row]][-1]
    for(m in mappings){
      results$Mapping[row] = paste0(results$Mapping[row], ", ",
                                   results$Mapping[results$MODULE==results$MODULE[row] &
                                                     results$Mapping==m], ": ",
                                   results$FDR[results$MODULE==results$MODULE[row] &
                                                     results$Mapping==m], " ")
    }
  }
}
results = results[!duplicated(results$MODULE),] # preserves FDR but not other information

# get gene length of modules
all_modules <- read.delim(paste0("./Mergeomics/", name, "_combined_module.txt"), stringsAsFactors = FALSE)
signif_modules = unique(results$MODULE)
geneNum = c()
for(mod in signif_modules){
  geneNum = append(geneNum, sum(all_modules$MODULE==mod))
}
results$ModuleGeneCount = geneNum
write.table(results, "MSEA_trimmed_compiled_results.txt", row.names = FALSE, quote = FALSE, sep = "\t")

selected_modules = data.frame(stringsAsFactors = FALSE)
for(mod in signif_modules){
  selected_modules = rbind(selected_modules, all_modules[all_modules$MODULE==mod, ])
}
colnames(selected_modules) = c("MODULE", "NODE")
write.table(selected_modules, "Significant_modules.txt", row.names = FALSE, quote = FALSE, sep = "\t")

## Key driver analysis

source("./utils/Mergeomics_YZ_2020_v2.R")
brainNet <- "./utils/brainNet_v1v3.txt"

job.kda <- list()
job.kda$label <- "KDA"
job.kda$folder <- "OUD"
job.kda$netfile <- brainNet
job.kda$modfile <- "Significant_modules.txt"

## "0" means we do not consider edge weights while 1 is opposite.
job.kda$edgefactor <- 1
## The searching depth for the KDA
job.kda$depth <- 1
## 0 means we do not consider the directions of the regulatory interactions
## while 1 is opposite.
job.kda$direction <- 0
job.kda$nperm <- 10000 # 10000 for formal analysis

moddata <- tool.read(job.kda$modfile)
mod.names <- unique(moddata$MODULE)
moddata <- moddata[which(!is.na(match(moddata$MODULE, mod.names))), ]
## save this to a temporary file and set its path as new job.kda$modfile:
tool.save(moddata, "subsetof.supersets.txt")
job.kda$modfile <- "subsetof.supersets.txt"

job.kda <- kda.configure(job.kda)
job.kda <- kda.start(job.kda)
job.kda <- kda.prepare(job.kda)
job.kda <- kda.analyze(job.kda)
job.kda <- kda.finish(job.kda)

job.kda <- kda2cytoscape(job.kda, ndrivers = 10)

# ============================================================
# Visualization: key-driver (KD) ranking prep for network figure
# ============================================================

### Modules

## THC

# Frequency of sig. KD occurence

kds.thc = read.csv("./WGCNA-RNA/Treatment/all-KDs.csv")
kds.thc = kds.thc[, -1]
kds.thc$Association = "THC"
kd.rank.thc = as.data.frame(table(kds.thc$NODE))
kd.rank.thc = kd.rank.thc[order(-kd.rank.thc$Freq), ]
write.csv(kd.rank.thc, "./WGCNA-RNA/Treatment/rank_all-KDs.csv")

# using top 10 KDs only

library(dplyr)
library(gtools)

top.kds.thc <- kds.thc %>% group_by(MODULE) %>% slice_min(FDR, n = 10)
top.kd.rank.thc = as.data.frame(table(top.kds.thc$NODE))
top.kd.rank.thc = top.kd.rank.thc[order(-top.kd.rank.thc$Freq), ]
write.csv(top.kd.rank.thc, "./WGCNA-RNA/Treatment/rank_top10-KDs.csv")

## Genotype

# Frequency of sig. KD occurence

kds.geno <- read.delim("./WGCNA-RNA/Genotype/Geno_Significant_modules/kda/KDA.results.txt")
kds.geno <- kds.geno[kds.geno$FDR < 0.05, ]
kds.geno$Association = "Geno"
write.csv(kds.geno, file = "./WGCNA-RNA/Genotype/Geno_Significant_modules/all-KDs.csv")

kd.rank.geno = as.data.frame(table(kds.geno$NODE))
kd.rank.geno = kd.rank.geno[order(-kd.rank.geno$Freq), ]
write.csv(kd.rank.geno, "./WGCNA-RNA/Genotype/Geno_Significant_modules/rank_all-KDs.csv")

# using top 10 KDs only

top.kds.geno <- kds.geno %>% group_by(MODULE) %>% slice_min(FDR, n = 10)
top.kd.rank.geno = as.data.frame(table(top.kds.geno$NODE))
top.kd.rank.geno = top.kd.rank.geno[order(-top.kd.rank.geno$Freq), ]
write.csv(top.kd.rank.geno, "./WGCNA-RNA/Genotype/Geno_Significant_modules/rank_top10-KDs.csv")

## THC + Genotype

# Frequency of sig. KD occurence

kds <- read.delim("./WGCNA-RNA/Treatment\ x\ Genotype/THCxHR_Significant_modules/kda/KDA.results.txt")
kds <- kds[kds$FDR < 0.05, ]
write.csv(kds, file = "./WGCNA-RNA/Treatment\ x\ Genotype/THCxHR_Significant_modules/all-KDs.csv")

kd.rank = as.data.frame(table(kds$NODE))
kd.rank = kd.rank[order(-kd.rank$Freq), ]
write.csv(kd.rank, "./WGCNA-RNA/Treatment\ x\ Genotype/THCxHR_Significant_modules/rank_all-KDs.csv")

# using top 10 KDs only

top.kds <- kds %>% group_by(MODULE) %>% slice_min(FDR, n = 10)
top.kd.rank = as.data.frame(table(top.kds$NODE))
top.kd.rank = top.kd.rank[order(-top.kd.rank$Freq), ]
write.csv(top.kd.rank, "./WGCNA-RNA/Treatment\ x\ Genotype/THCxHR_Significant_modules/rank_top10-KDs.csv")

# rank using unique(c(kds.thc$MODULE, kds.geno$MODULE))

shared.mod = intersect(kds.thc$MODULE, kds.geno$MODULE)
all.kds = rbind(kds.thc, subset(kds.geno, !(MODULE %in% shared.mod)))

all.kd.rank = as.data.frame(table(all.kds$NODE))
all.kd.rank = all.kd.rank[order(-all.kd.rank$Freq), ]

all.top.kds <- all.kds %>% group_by(MODULE) %>% slice_min(FDR, n = 10)
all.top.kd.rank = as.data.frame(table(all.top.kds$NODE))
all.top.kd.rank = all.top.kd.rank[order(-all.top.kd.rank$Freq), ]

### DEG's KDs

## THC

# Frequency of sig. KD occurence

deg.kds.thc  <- read.delim("./Merged/DEGs/Seurat-DEGs_no-ribo/THC-DEG_KDA/kda/KDA.results.txt")
deg.kds.thc <- deg.kds.thc[deg.kds.thc$FDR < 0.05, ]
deg.kds.thc$Association = "THC"
write.csv(deg.kds.thc, file = "./Merged/DEGs/Seurat-DEGs_no-ribo/THC-DEG_KDA/all-KDs.csv")

deg.kd.rank.thc = as.data.frame(table(deg.kds.thc$NODE))
deg.kd.rank.thc = deg.kd.rank.thc[order(-deg.kd.rank.thc$Freq), ]
write.csv(deg.kd.rank.thc, "./Merged/DEGs/Seurat-DEGs_no-ribo/THC-DEG_KDA/rank_all-KDs.csv")

# using top 10 KDs only

deg.top.kds.thc <- deg.kds.thc %>% group_by(MODULE) %>% slice_min(FDR, n = 10)
deg.top.kd.rank.thc = as.data.frame(table(deg.top.kds.thc$NODE))
deg.top.kd.rank.thc = deg.top.kd.rank.thc[order(-deg.top.kd.rank.thc$Freq), ]
write.csv(deg.top.kd.rank.thc, "./Merged/DEGs/Seurat-DEGs_no-ribo/THC-DEG_KDA/rank_top10-KDs.csv")

## Geno

# Frequency of sig. KD occurence

deg.kds.geno  <- read.delim("./Merged/DEGs/Seurat-DEGs_no-ribo/Geno-DEG_KDA/kda/KDA.results.txt")
deg.kds.geno <- deg.kds.geno[deg.kds.geno$FDR < 0.05, ]
deg.kds.geno$Association = "Geno"
write.csv(deg.kds.thc, file = "./Merged/DEGs/Seurat-DEGs_no-ribo/Geno-DEG_KDA/all-KDs.csv")

deg.kd.rank.geno = as.data.frame(table(deg.kds.geno$NODE))
deg.kd.rank.geno = deg.kd.rank.geno[order(-deg.kd.rank.geno$Freq), ]
write.csv(deg.kd.rank.geno, "./Merged/DEGs/Seurat-DEGs_no-ribo/Geno-DEG_KDA/rank_all-KDs.csv")

### All of them

all.kds = rbind(all.kds, deg.kds.thc, deg.kds.geno)

all.kd.rank = as.data.frame(table(all.kds$NODE))
all.kd.rank = all.kd.rank[order(-all.kd.rank$Freq), ]
all.kd.rank$Association = apply(all.kd.rank, 1, function(x) paste0(all.kds$Association[all.kds$NODE == x[1]], collapse = ", "))
all.kd.rank$source = apply(all.kd.rank, 1, function(x) paste0(all.kds$MODULE[all.kds$NODE == x[1]], collapse = ", "))

inter.kds.rank = all.kd.rank %>% filter(grepl('THC, Geno|Geno, THC', Association))
inter.kds = all.kds %>% filter(NODE %in% inter.kds.rank$Var1) # 289

all.top.kds <- all.kds %>% group_by(MODULE) %>% slice_min(FDR, n = 10)
all.top.kd.rank = as.data.frame(table(all.top.kds$NODE))
all.top.kd.rank = all.top.kd.rank[order(-all.top.kd.rank$Freq), ]
all.top.kd.rank$Association = apply(all.top.kd.rank, 1, function(x) paste0(all.top.kds$Association[all.top.kds$NODE == x[1]], collapse = ", "))
all.top.kd.rank$Source = apply(all.top.kd.rank, 1, function(x) paste0(all.top.kds$MODULE[all.top.kds$NODE == x[1]], collapse = ", "))

inter.top.kds.rank = all.top.kd.rank %>% filter(grepl('THC, Geno|Geno, THC', Association)) # 25
inter.top.kds = all.top.kds %>% filter(NODE %in% inter.top.kds.rank$Var1)

## THC

# Frequency of sig. KD occurence

all.kds.thc <- all.kds %>% filter(Association == "THC")

all.kd.rank.thc = as.data.frame(table(all.kds.thc$NODE))
all.kd.rank.thc = all.kd.rank.thc[order(-all.kd.rank.thc$Freq), ]

# using top 10 KDs only

all.top.kds.thc <- all.kds.thc %>% group_by(MODULE) %>% slice_min(FDR, n = 10)
all.top.kd.rank.thc = as.data.frame(table(all.top.kds.thc$NODE))
all.top.kd.rank.thc = all.top.kd.rank.thc[order(-all.top.kd.rank.thc$Freq), ]

## Geno

# Frequency of sig. KD occurence

all.kds.geno <- all.kds %>% filter(Association == "Geno")

all.kd.rank.geno = as.data.frame(table(all.kds.geno$NODE))
all.kd.rank.geno = all.kd.rank.geno[order(-all.kd.rank.geno$Freq), ]
all.kd.rank.geno$source = apply(all.kd.rank.geno, 1, function(x) paste0(all.kds.geno$MODULE[all.kds.geno$NODE == x[1]], collapse = ", "))

# using top 10 KDs only

all.top.kds.geno <- all.kds.geno %>% group_by(MODULE) %>% slice_min(FDR, n = 10)
all.top.kd.rank.geno = as.data.frame(table(all.top.kds.geno$NODE))
all.top.kd.rank.geno = all.top.kd.rank.geno[order(-all.top.kd.rank.geno$Freq), ]
all.top.kd.rank.geno$source = apply(all.top.kd.rank.geno, 1, function(x) paste0(all.kds.geno$MODULE[all.kds.geno$NODE == x[1]], collapse = ", "))

all.kds$Association[which(all.kds$MODULE %in% shared.mod)] = "THC, Geno"

# ============================================================
# Visualization: mergeomics MSEA summary heatmap across traits
# ============================================================

library(pheatmap)
library(RColorBrewer)
library(stringr)

cor.sig <- read.csv("./WGCNA-RNA/all-sig_moduleTraitCor.csv")
mod.sig <- read.csv("./WGCNA-RNA/all-sig_moduleTraitPvalue.csv")
cor.sig <- cor.sig[, -1]; mod.sig <- mod.sig[, -1]

sig.mod <- mod.sig$X # 43 modules

deg <- read.delim("./Merged/Seurat-DEGs_no-ribo/DEG-HGNC_treatment.txt")
sig.mod <- c(sig.mod, unique(deg$MODULE))

deg <- read.delim("./Merged/Seurat-DEGs_no-ribo/DEG-HGNC_genotype.txt")
sig.mod <- c(sig.mod, unique(deg$MODULE))

cud <- read.delim ("./Mergeomics/CUD/NAc-module-DEG/MSEA_trimmed_compiled_results.txt")
cud <- cud[which(cud$MODULE %in% sig.mod), ]

oud <- read.delim ("./Mergeomics/OUD/NAc-module-DEG/MSEA_trimmed_compiled_results.txt")
oud <- oud[which(oud$MODULE %in% sig.mod), ]
oud$Mapping <- paste0("oud-", oud$Mapping)

merg.sum.sigmods <- rbind(cud, oud)

scz <- read.delim ("./Mergeomics/SCZ/NAc-module-DEG/MSEA_trimmed_compiled_results.txt")
scz <- scz[which(scz$MODULE %in% sig.mod), ]
scz$Mapping <- paste0("scz-", scz$Mapping)

merg.sum.sigmods <- rbind(merg.sum.sigmods, scz)

ga <- read.delim ("./Mergeomics/General-addiction/NAc-module-DEG/MSEA_trimmed_compiled_results.txt")
ga <- ga[which(ga$MODULE %in% sig.mod), ]
ga$Mapping <- paste0("ga-", ga$Mapping)

merg.sum.sigmods <- rbind(merg.sum.sigmods, ga)

agesmk <- read.delim ("./Mergeomics/TAUD_EUR/AgeSmk/NAc-module-DEG/MSEA_trimmed_compiled_results.txt")
agesmk <- agesmk[which(agesmk$MODULE %in% sig.mod), ]
agesmk$Mapping <- paste0("agesmk-", agesmk$Mapping)

merg.sum.sigmods <- rbind(merg.sum.sigmods, agesmk)

cigday <- read.delim ("./Mergeomics/TAUD_EUR/CigDay/NAc-module-DEG/MSEA_trimmed_compiled_results.txt")
cigday <- cigday[which(cigday$MODULE %in% sig.mod), ]
cigday$Mapping <- paste0("cigday-", cigday$Mapping)

merg.sum.sigmods <- rbind(merg.sum.sigmods, cigday)

DrnkWk <- read.delim ("./Mergeomics/TAUD_EUR/DrnkWk/NAc-module-DEG/MSEA_trimmed_compiled_results.txt")
DrnkWk <- DrnkWk[which(DrnkWk$MODULE %in% sig.mod), ]
DrnkWk$Mapping <- paste0("DrnkWk-", DrnkWk$Mapping)

merg.sum.sigmods <- rbind(merg.sum.sigmods, DrnkWk)

SmkCes <- read.delim ("./Mergeomics/TAUD_EUR/SmkCes/NAc-module-DEG/MSEA_trimmed_compiled_results.txt")
SmkCes <- SmkCes[which(SmkCes$MODULE %in% sig.mod), ]
SmkCes$Mapping <- paste0("SmkCes-", SmkCes$Mapping)

merg.sum.sigmods <- rbind(merg.sum.sigmods, SmkCes)

SmkInit <- read.delim ("./Mergeomics/TAUD_EUR/SmkInit/NAc-module-DEG/MSEA_trimmed_compiled_results.txt")
SmkInit <- SmkInit[which(SmkInit$MODULE %in% sig.mod), ]
SmkInit$Mapping <- paste0("SmkInit-", SmkInit$Mapping)

merg.sum.sigmods <- rbind(merg.sum.sigmods, SmkInit)

ADHD <- read.delim ("./Mergeomics/ADHD/NAc-module-DEG/MSEA_trimmed_compiled_results.txt")
ADHD <- ADHD[which(ADHD$MODULE %in% sig.mod), ]
ADHD$Mapping <- paste0("ADHD-", ADHD$Mapping)

merg.sum.sigmods <- rbind(merg.sum.sigmods, ADHD)

Externalizing <- read.delim ("./Mergeomics/Externalizing/NAc-module-DEG/MSEA_trimmed_compiled_results.txt")
Externalizing <- Externalizing[which(Externalizing$MODULE %in% sig.mod), ]
Externalizing$Mapping <- paste0("Externalizing-", Externalizing$Mapping)

merg.sum.sigmods <- rbind(merg.sum.sigmods, Externalizing)

bmd <- read.delim ("./Mergeomics/BMD/NAc-module-DEG/MSEA_trimmed_compiled_results.txt")
bmd <- bmd[which(bmd$MODULE %in% sig.mod), ]
bmd$Mapping <- paste0("BMD-", bmd$Mapping)

merg.sum.sigmods <- rbind(merg.sum.sigmods, bmd)

lbm <- read.delim ("./Mergeomics/LBM/NAc-module-DEG/MSEA_trimmed_compiled_results.txt")
lbm <- lbm[which(lbm$MODULE %in% sig.mod), ]
lbm$Mapping <- paste0("LBM-", lbm$Mapping)

merg.sum.sigmods <- rbind(merg.sum.sigmods, lbm)

# mannually put together and organize

# Significant modules
sum.sigmod <- read.csv("./Mergeomics/summary_sig.mod.csv")
rownames(sum.sigmod) <- sum.sigmod[, 1]
sum.sigmod <- sum.sigmod[, -1]
sum.sigmod <- t(sum.sigmod)

test_labels <- sum.sigmod

test_labels <- as.matrix(test_labels)
test_labels[test_labels < 0.001] <- "***"
test_labels[test_labels < 0.01 & test_labels >= 0.001] <- "**"
test_labels[test_labels < 0.05 & test_labels >= 0.01] <- "*"
test_labels[test_labels >= 0.05] <- ""

sum.sigmod <- -log10(sum.sigmod)

pheatmap(sum.sigmod, scale = "none", trace = "none", col = colorRampPalette(brewer.pal(9, "Reds"))(20), cluster_cols = FALSE, cluster_rows = FALSE, border_color = NA, display_numbers = test_labels, fontsize_number = 30) # 12.5 x 8.33

# ============================================================
# Visualization: Bayesian network node/edge styling for Cytoscape
# (originally Scripts/20250220_BN.txt)
#
# Matches manuscript Fig. 4E (Zuo et al., bioRxiv 2025.04.03.646846):
# "Bayesian network representation of the NB-1_blue module. Diamond
# nodes denote key driver genes; hexagons denote hub genes; blue
# nodes represent genes within the NB-1_blue module; gray nodes
# represent genes present in the Bayesian network but not assigned
# to these categories." The final network image itself is rendered
# in Cytoscape (v3.8.2) from the node/edge tables produced below;
# this script only computes node color/shape/size attributes.
#
# Note: this script was saved using the module's pre-rename label,
# GABA-2 (blue module) -- GABA-2 was renamed to NB-1 later in the
# pipeline (see 20250324_rename.txt), so GABA-2_blue = NB-1_blue.
# The hub gene list and key-driver network below were computed for
# that module specifically; re-run against current file names/paths
# if the underlying KDA/WGCNA outputs have since moved.
# ============================================================

# Visualize SATB1, NR3C2, CNR1 (example key drivers) in the NB-1_blue
# (formerly GABA-2_blue) Bayesian network

nodes <- read.delim("./WGCNA-RNA/Treatment\ x\ Genotype/THCxHR_Significant_modules-nKD50/cytoscape/kda2cytoscape.nodes.txt", stringsAsFactors = FALSE)

mod <- read.delim("./WGCNA-RNA/HGNC_outputs/GABA-2_WGCNA.txt", stringsAsFactors = FALSE)
mod <- mod[which(mod$MODULE == "blue"), ]
blue = mod

mod <- read.delim("./Merged/Seurat-DEGs_no-ribo/DEG-HGNC_treatment.txt", stringsAsFactors = FALSE)
mod <- mod[which(mod$MODULE == "GABA-2_DEG-Treatment"), ]
deg = mod
dups <- nodes$NODE[which(nodes$NODE %in% blue$GENE & nodes$NODE %in% deg$NODE)]

load("./Merged/DEGs/20230803_DEGs/DEG.RData")
deg = sig.thc.markers[["GABA-2 (Tiam2+)"]]
deg.up = deg$Gene[which(deg$avg_log2FC > 0)]
deg.down = deg$Gene[which(deg$avg_log2FC < 0)]
deg.reln = sig.reln.markers[["GABA-2 (Tiam2+)"]]$Gene

library(dplyr)
mouse_human_genes = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt", sep = "\t")
# write.csv(mouse_human_genes, "./utils/mouse-human-genes.csv")

convert_mouse_to_human <- function(gene_list){

  require("dplyr")
  output = data.frame(Input = gene_list, Output = rep("", length(gene_list)))

  for(gene in gene_list){
    class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name == "mouse, laboratory"))[['DB.Class.Key']]
    if(!identical(class_key, integer(0))){
      human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name == "human"))[, "Symbol"]
      for(human_gene in human_genes){
        output$Output[which(output$Input == gene)] = human_gene
      }
    }
  }
  no.gene = which(output$Output == "")
  output$Output[no.gene] = output$Input[no.gene]
  return(output)
}

convert_human_to_mouse <- function(gene_list){

  require("dplyr")
  mouse_human_genes = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt", sep = "\t")
  output = data.frame(Input = gene_list, Output = rep("", length(gene_list)))

  for(gene in gene_list){
    class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name == "human"))[['DB.Class.Key']]
    if(!identical(class_key, integer(0))){
      mouse_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name == "mouse, laboratory"))[, "Symbol"]
      for(mouse_gene in mouse_genes){
        output$Output[which(output$Input == gene)] = mouse_gene
      }
    }
  }
  no.gene = which(output$Output == "")
  output$Output[no.gene] = output$Input[no.gene]
  return(output)
}

deg.up = convert_mouse_to_human(deg.up)
deg.down = convert_mouse_to_human(deg.down)
deg.reln = convert_mouse_to_human(deg.reln)

# node color: module membership (blue), up-regulated THC DEGs, down-regulated THC DEGs
# https://projects.susielu.com/viz-palette?colors=[%22#93a3fa%22,%22#ffa8a0%22,%22#a8f0a1%22,%22#ededed%22]&backgroundColor=%22white%22&fontColor=%22black%22&mode=%22normal%22

nodes$COLOR <- rep("#ededed", length(nodes$NODE))
nodes$COLOR[which(nodes$NODE %in% blue$GENE)] <- "#93A3FA"
nodes$COLOR[which(nodes$NODE %in% deg.up$Output)] <- "#FFA8A0"
nodes$COLOR[which(nodes$NODE %in% deg.down$Output)] <- "#A8F0A1"

nodes$Mouse <- nodes$LABEL
nodes$Mouse <- convert_human_to_mouse(nodes$Mouse)$Output

# node shape: hexagon = WGCNA hub genes, diamond = key driver (Bayesian network) genes
nodes$SHAPE <- rep("Ellipse", length(nodes$NODE))
hubs <- c("Nfib", "Tiam2", "Nol4", "Zbtb20", "Dlx6os1", "Ctnna2", "Enox1", "Nfia", "Plcl1", "Zfp704", "Nav3", "Ccnd2", "Msi2", "Igfbpl1", "Mpped2", "Plxna4", "CT025619.1", "Stmn2", "Spon1", "Nrxn3", "Pcdh9", "Tcf4", "Gm26936", "Ccdc88a", "Gm3764")
nodes$SHAPE[which(nodes$Mouse %in% hubs)] <- "Hexagon"

edges <- read.delim("./WGCNA-RNA/Treatment\ x\ Genotype/THCxHR_Significant_modules-nKD50/cytoscape/kda2cytoscape.edges.txt", stringsAsFactors = FALSE)
bn.hubs <- unique(edges$TAIL)
nodes$SHAPE[which(nodes$NODE %in% bn.hubs)] <- "Diamond"

# node size and membership flags
nodes$SIZE[which(nodes$SIZE == 50)] <- 10
nodes$MEMBERSHIP <- rep("NONMEMBER", length(nodes$NODE))
nodes$MEMBERSHIP[which(nodes$NODE %in% blue$GENE | nodes$NODE %in% c(deg.up$Output, deg.down$Output))] <- "MEMBER"
nodes$SIZE[which(nodes$COLOR != "#ededed" & nodes$SIZE == 10)] <- 100
nodes$SIZE[which(nodes$COLOR == "#ededed" & nodes$SIZE != 200)] <- 50
nodes$SIZE[which(nodes$NODE %in% setdiff(bn.hubs, c("SATB1", "NR3C2", "CNR1")))] <- 150
nodes$FONT[which(nodes$NODE %in% setdiff(bn.hubs, c("SATB1", "NR3C2", "CNR1")))] <- 50

nodes$width <- rep(1, length(nodes$NODE))
nodes$width[which(nodes$Mouse %in% hubs | nodes$NODE %in% bn.hubs)] <- 40
nodes$width[which(nodes$width == 40 | nodes$SIZE == 100)] <- 20
nodes$width[which(nodes$width == 20 & nodes$SIZE == 200)] <- 40
nodes$width[which(nodes$width == 40 & nodes$SHAPE == "Ellipse")] <- 20

# border color: flag genes also differentially expressed by Reln genotype
nodes$border <- rep("#ffffff", length(nodes$NODE))
nodes$border[which(nodes$NODE %in% deg.reln$Output)] <- "#f5ac40"

write.table(nodes, "./WGCNA-RNA/gaba2-blue_nodes.txt", sep = "\t",quote = FALSE,row.names = FALSE, col.names = TRUE)

edges <- edges[which(edges$TAIL %in% c("SATB1", "NR3C2", "CNR1")), ]
write.table(edges, "./WGCNA-RNA/gaba2-blue_edges.txt", sep = "\t", quote = FALSE,row.names = FALSE, col.names = TRUE)

# Venn diagram: overlap between NB-1_blue module genes and THC DEGs in this
# cell type (a module-vs-DEG sanity check; not the Fig. 3E Factor-3-vs-DEG
# Venn, which is part of the colleague-authored Tensor-Cell2Cell
# communication analysis and out of scope for this repo)
library(VennDiagram)
library(RColorBrewer)
vd <- VennDiagram::venn.diagram(list(Module = blue$GENE, DEG = deg$NODE), filename = NULL)
grid::grid.draw(vd)
