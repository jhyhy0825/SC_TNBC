#! /usr/bin/env Rscript
  
# Updated by Hye-Yeon Ju on Apr 22, 2024.

##### Code Description #####
# This code performs high dimensional weighted gene co-expression network analysis using hdWGCNA.

##### Code Example #####
# Rscript
# /home/jhy/single_cell_TNBC/code/for_github/4_High_dimensional_WGCNA/hdWGCNA.R

library(Seurat)

### Import RDS file ###
data <- readRDS('/home/jhy/single_cell_TNBC/result/Seurat/integration2/result2/Integration_mt10_feat4000_dim10_res0.4_with_case_control.rds')

### Specify clusters based on the cell type annotation result ##
B <- subset(data, seurat_clusters %in% c(5,12,15))
T <- subset(data, seurat_clusters %in% c(0,1,3,6,10,11,14,16))
NK <- subset(data, seurat_clusters %in% c(2))
Mye <- subset(data, seurat_clusters %in% c(4,7,8,9,13,17))
raw_seurat_list <- c(T,B,NK,Mye)
merged_seurat <- merge(x = raw_seurat_list[[1]],y = raw_seurat_list[2:length(raw_seurat_list)],merge.data = TRUE)
merged_seurat$umap <- data$umap
merged_seurat$pca <- data$pca
merged_seurat$harmony <- data$harmony

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

###Tcell불러온다음에###
##level확인 중요함 순서가 0,1,2 아닐수도 있음##
data <- readRDS('/home/jhy/single_cell_TNBC/result/Seurat/integration2/result2/Integration_mt10_feat4000_dim10_res0.4_tcell.rds')
new.cluster.ids <- c('Central memory T cell', 'Terminally differentiated T cell','Effector memory T cell','Terminally differentiated T cell','Effector memory T cell','Terminally differentiated T cell','Terminally differentiated T cell','Terminally differentiated T cell','Terminally differentiated T cell','Effector memory T cell', 'Terminally differentiated T cell')
names(new.cluster.ids) <- levels(data)
data <- RenameIdents(data, new.cluster.ids)
data$cell_type <- Idents(data)
###

##
ori <- readRDS('/home/jhy/single_cell_TNBC/result/Seurat/integration2/result2/Integration_mt10_feat4000_dim10_res0.4_with_case_control.rds')
data <- ori
Idents(data) <- 'seurat_clusters'
new.cluster.ids <- c('NK cell','T cell','T cell','T cell','B cell','T cell','Myeloid cell','B cell','Myeloid cell','B cell','T cell','T cell','Myeloid cell','Myeloid cell','Myeloid cell','Myeloid cell','T cell','T cell')
names(new.cluster.ids) <- levels(data)
data <- RenameIdents(data, new.cluster.ids)
data$cell_type <- Idents(data)
##

seurat_obj <- data
seurat_obj <- SetupForWGCNA(seurat_obj,gene_select = "fraction",fraction = 0.05,wgcna_name = "tutorial")
seurat_obj <- MetacellsByGroups(seurat_obj = seurat_obj,group.by = c("cell_type"),reduction = 'harmony',k = 25,max_shared = 10,ident.group = 'cell_type')
### groupby에celltype?
seurat_obj <- NormalizeMetacells(seurat_obj)
#seurat_obj <- SetDatExpr(seurat_obj,group_name = "T cell",group.by='cell_type',assay = 'RNA',slot = 'data')
seurat_obj <- SetDatExpr(seurat_obj,group_name = "T cell",group.by='cell_type',assay = 'SCT',slot = 'data')
seurat_obj <- TestSoftPowers(seurat_obj,networkType = 'signed')
pdf('/home/jhy/single_cell_TNBC/result/Seurat/integration2/result2/test/hdWGCNA_softpower.pdf')
plot_list <- PlotSoftPowers(seurat_obj)
wrap_plots(plot_list, ncol=2)
dev.off()
power_table <- GetPowerTable(seurat_obj)
head(power_table)

seurat_obj <- ConstructNetwork(seurat_obj, tom_name = 'T cell')
pdf('/home/jhy/single_cell_TNBC/result/Seurat/integration2/result2/test/hdWGCNA_dendrogrma.pdf')
PlotDendrogram(seurat_obj, main='INH hdWGCNA Dendrogram')
dev.off()

seurat_obj <- ModuleEigengenes(seurat_obj)
hMEs <- GetMEs(seurat_obj)
MEs <- GetMEs(seurat_obj, harmonized=FALSE)
seurat_obj <- ModuleConnectivity(seurat_obj,group.by = 'cell_type', group_name = 'T cell')
seurat_obj <- ResetModuleNames(seurat_obj,new_name = "Tcell-M")
pdf('/home/jhy/single_cell_TNBC/result/Seurat/integration2/result2/test/hdWGCNA_module.pdf')
p <- PlotKMEs(seurat_obj, ncol=5)
p
dev.off()

### Visualization ###
pdf('/home/jhy/single_cell_TNBC/result/Seurat/integration2/result2/test/hdWGCNA_with_UMAP.pdf')
plot_list <- ModuleFeaturePlot(seurat_obj,features='hMEs',order=TRUE)
wrap_plots(plot_list, ncol=6)
dev.off()

pdf('/home/jhy/single_cell_TNBC/result/Seurat/integration2/result2/test/hdWGCNA_correlation.pdf')
ModuleCorrelogram(seurat_obj)
dev.off()

