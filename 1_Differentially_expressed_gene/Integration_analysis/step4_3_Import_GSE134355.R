#! /usr/bin/env Rscript

# Updated by Hye-Yeon Ju on Nov 24, 2023.

##### Code Description #####
# This code imports public data (GSE134355, healthy donors)

##### Code Example #####
# Rscript
## /home/jhy/single_cell_TNBC/code/integration/For_GSE134355.R

library(dplyr)
library(Seurat)
library(patchwork)

outputdir <- '/home/jhy/single_cell_TNBC/result/Seurat/integration2/'

data1 <- read.table(gzfile('/home/jhy/single_cell_TNBC/result/Seurat/integration2/GSE134355/new/AdultPeripheralBlood1.rmbatchdge.txt.gz'),sep=' ',header=T,row.names=1)
data2 <- read.table(gzfile('/home/jhy/single_cell_TNBC/result/Seurat/integration2/GSE134355/new/AdultPeripheralBlood2.rmbatchdge.txt.gz'),sep=' ',header=T,row.names=1)
data3 <- read.table(gzfile('/home/jhy/single_cell_TNBC/result/Seurat/integration2/GSE134355/new/AdultPeripheralBlood3.rmbatchdge.txt.gz'),sep=' ',header=T,row.names=1)
data4 <- read.table(gzfile('/home/jhy/single_cell_TNBC/result/Seurat/integration2/GSE134355/new/AdultPeripheralBlood4.rmbatchdge.txt.gz'),sep=' ',header=T,row.names=1)

data1 <- CreateSeuratObject(counts = data1, project = 'Normal1', min.cells = 3, min.features = 200)
data1
data2 <- CreateSeuratObject(counts = data2, project = 'Normal2', min.cells = 3, min.features = 200)
data2
data3 <- CreateSeuratObject(counts = data3, project = 'Normal3', min.cells = 3, min.features = 200)
data3
data4 <- CreateSeuratObject(counts = data4, project = 'Normal4', min.cells = 3, min.features = 200)
data4

pbmc <- merge(data1, y = c(data2,data3,data4), add.cell.ids = c('Normal1','Normal2','Normal3','Normal4'), project = 'GSE134355')
pbmc


pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")
# Visualize_QC_metrics
pdf(file=paste0(outputdir,"feature_vlnplot_GSE134355.pdf"))
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
# Visualize_feature-feature_relationships
pdf(file=paste0(outputdir,"feature_scatter_GSE134355.pdf"),width=13,height=10)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()

saveRDS(pbmc, file = paste0(outputdir,"/preprocessing_GSE134355.rds"))

