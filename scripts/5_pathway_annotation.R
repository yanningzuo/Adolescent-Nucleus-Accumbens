# Pathway analysis
library(enrichR)
dbs <- c("KEGG_2021_Human", "GO_Biological_Process_2021", "BioCarta_2016", "Reactome_2022")
reln.deg.go.up <- list()
for (clu in levels(Idents(data))){
    temp <- sig.reln.markers[[clu]]
    temp <- temp[which(temp$avg_log2FC > 0), ]
    genes <- rownames(temp)
    go.genes <- enrichr(genes, dbs)
    if (is.null(go.genes) == TRUE){
      reln.deg.go.up[[clu]] <- 0
    }else{
      list.go.genes = do.call("rbind",go.genes)
      list.go.genes$padj <- p.adjust(list.go.genes$P.value, method = "BH")
      reln.deg.go.up[[clu]] <- list.go.genes[order(list.go.genes$padj), ]
    }
}
reln.deg.go.down <- list()
for (clu in levels(Idents(data))){
    temp <- sig.reln.markers[[clu]]
    temp <- temp[which(temp$avg_log2FC < 0), ]
    genes <- rownames(temp)
    go.genes <- enrichr(genes, dbs)
    if (is.null(go.genes) == TRUE){
      reln.deg.go.down[[clu]] <- 0
    }else{
      list.go.genes = do.call("rbind", go.genes)
      list.go.genes$padj <- p.adjust(list.go.genes$P.value, method = "BH")
      reln.deg.go.down[[clu]] <- list.go.genes[order(list.go.genes$padj), ]
    }
}
thc.deg.go.up <- list()
for (clu in levels(Idents(data))){
    temp <- sig.thc.markers[[clu]]
    temp <- temp[which(temp$avg_log2FC > 0), ]
    genes <- rownames(temp)
    go.genes <- enrichr(genes, dbs)
    if (is.null(go.genes) == TRUE){
      thc.deg.go.up[[clu]] <- 0
    }else{
      list.go.genes = do.call("rbind", go.genes)
      list.go.genes$padj <- p.adjust(list.go.genes$P.value, method = "BH")
      thc.deg.go.up[[clu]] <- list.go.genes[order(list.go.genes$padj), ]
    }
}
thc.deg.go.down <- list()
for (clu in levels(Idents(data))){
    temp <- sig.thc.markers[[clu]]
    temp <- temp[which(temp$avg_log2FC < 0), ]
    genes <- rownames(temp)
    go.genes <- enrichr(genes, dbs)
    if (is.null(go.genes) == TRUE){
      thc.deg.go.down[[clu]] <- 0
    }else{
      list.go.genes = do.call("rbind", go.genes)
      list.go.genes$padj <- p.adjust(list.go.genes$P.value, method = "BH")
      thc.deg.go.down[[clu]] <- list.go.genes[order(list.go.genes$padj), ]
    }
}
save.image("DEG.RData")
