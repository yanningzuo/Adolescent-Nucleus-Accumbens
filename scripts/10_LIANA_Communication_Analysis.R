library(tidyverse)
library(magrittr)
library(liana)
library(OmnipathR)
library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)
set.seed(12)
attach('20230107_SCT_noGlut-61911.RData')
data <- data
detach('file:20230107_SCT_noGlut-61911.RData')
DefaultAssay(data) <- "SCT"
Idents(object = data) <- data@meta.data[["detailed.celltype"]]
# convert LIANA's Consensus resource to murine symbols
op_resource <- select_resource("Consensus")[[1]]
# generate orthologous resource
ortholog_resource <- generate_homologs(op_resource = op_resource,
                                       target_organism = 10090) # mouse
# run LIANA with the orthologous resource
liana_res <- liana_wrap(data, resource = 'custom', external_resource = ortholog_resource)
# aggregate
liana_res <- liana_res %>%
    liana_aggregate()
# change name
liana_res$source[which(liana_res$source == "GABA-1 (Cdh18+)")] = "In-1"
liana_res$source[which(liana_res$source == "GABA-2 (Tiam2+)")] = "NB-1"
liana_res$source[which(liana_res$source == "GABA-3 (St18+ Kit+)")] = "In-2"
liana_res$source[which(liana_res$source == "GABA-4 (Cdh9+)")] = "In-3"
liana_res$source[which(liana_res$source == "GABA-5 (Vip+ Cck+)")] = "In-4"
liana_res$source[which(liana_res$source == "GABA-6 (Nos1+ Npy+)")] = "In-5"
liana_res$source[which(liana_res$source == "GABA-7 (Calb2+)")] = "NB-2"
liana_res$source[which(liana_res$source == "GABA-8 (Chat+)")] = "In-6"
liana_res$target[which(liana_res$target == "GABA-1 (Cdh18+)")] = "In-1"
liana_res$target[which(liana_res$target == "GABA-2 (Tiam2+)")] = "NB-1"
liana_res$target[which(liana_res$target == "GABA-3 (St18+ Kit+)")] = "In-2"
liana_res$target[which(liana_res$target == "GABA-4 (Cdh9+)")] = "In-3"
liana_res$target[which(liana_res$target == "GABA-5 (Vip+ Cck+)")] = "In-4"
liana_res$target[which(liana_res$target == "GABA-6 (Nos1+ Npy+)")] = "In-5"
liana_res$target[which(liana_res$target == "GABA-7 (Calb2+)")] = "NB-2"
liana_res$target[which(liana_res$target == "GABA-8 (Chat+)")] = "In-6"
liana_trunc <- liana_res %>%
  filter(source %in% c("In-4")) %>%
  filter(target == "NB-1") %>%
  filter(aggregate_rank <= 0.01)
