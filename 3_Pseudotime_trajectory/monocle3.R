#! /usr/bin/env Rscript
  
# Updated by Hye-Yeon Ju on Apr 22, 2024.

##### Code Description #####
# This code performs pseudotime trajectory analysis using monocle3.

##### Code Example #####
# Rscript
# /home/jhy/single_cell_TNBC/code/for_github/3_Pseudotime_trajectory/monocle3.R

library(Seurat)
library(tidyverse)
library(dplyr)
library(monocle3)
library(SeuratWrappers)
library(ggplot2)
library(ggridges)

### Import RDS file ###
data <- readRDS('/home/jhy/single_cell_TNBC/result/Seurat/integration2/result2/Integration_mt10_feat4000_dim10_res0.4.rds')
cds <- as.cell_data_set(data)
fData(cds)$gene_short_name <- rownames(fData(cds))
recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
recreate.partitions

### Pseudotime analysis ###
cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions
list.cluster <- data@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- data@reductions$umap@cell.embeddings
cds <- learn_graph(cds, use_partition = F)

### Visualization ###
pdf('/home/jhy/single_cell_TNBC/result/Seurat/integration2/result2/monocle3_tcell.pdf')
plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F,label_branch_points = T, label_roots = T, label_leaves = F,group_label_size = 5)
dev.off()

### Visualiation with new root ###
cds <- order_cells(cds, reduction_method = "UMAP", root_cells = colnames(cds[, clusters(cds) == 5]))
pdf('/home/jhy/single_cell_TNBC/result/Seurat/integration2/result2/monocle3_tcell_root.pdf')
plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = T,label_branch_points = T, label_roots = F, label_leaves = F)
dev.off()

