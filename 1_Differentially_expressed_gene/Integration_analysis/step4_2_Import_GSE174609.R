#! /usr/bin/env Rscript

# Updated by Hye-Yeon Ju on Nov 23, 2023.

##### Code Description #####
# This code imports public data (GSE174609, healthy donors)

##### Code Example #####
# Rscript
# /home/jhy/single_cell_TNBC/code/for_github/1_Differentially_expressed_gene/Integration_analysis/step4_2_Import_GSE174609.R

library(dplyr)
library(Seurat)
library(patchwork)

outputdir <- '/home/jhy/single_cell_TNBC/result/Seurat/integration2/'

#a <- read.table(gzfile(paste0(inputdir,'GSM4416534_PT-3232.csv.gz')),as.is=T,sep=',',header=T,row.names=1)

data1 <- read.table(gzfile('/home/jhy/single_cell_TNBC/result/Seurat/integration2/GSE174609/GSM5320459_Ctrl1_count_matrix.csv.gz'),sep=',',header=T,row.names=1)
data2 <- read.table(gzfile('/home/jhy/single_cell_TNBC/result/Seurat/integration2/GSE174609/GSM5320460_Ctrl2_count_matrix.csv.gz'),sep=',',header=T,row.names=1)
data3 <- read.table(gzfile('/home/jhy/single_cell_TNBC/result/Seurat/integration2/GSE174609/GSM5320461_Ctrl3_count_matrix.csv.gz'),sep=',',header=T,row.names=1)
data4 <- read.table(gzfile('/home/jhy/single_cell_TNBC/result/Seurat/integration2/GSE174609/GSM5320462_Ctrl4_count_matrix.csv.gz'),sep=',',header=T,row.names=1)

data1 <- CreateSeuratObject(counts = data1, project = 'Ctrl1', min.cells = 3, min.features = 200)
data1
data2 <- CreateSeuratObject(counts = data2, project = 'Ctrl2', min.cells = 3, min.features = 200)
data2
data3 <- CreateSeuratObject(counts = data3, project = 'Ctrl3', min.cells = 3, min.features = 200)
data3
data4 <- CreateSeuratObject(counts = data4, project = 'Ctrl4', min.cells = 3, min.features = 200)
data4

pbmc <- merge(data1, y = c(data2,data3,data4), add.cell.ids = c('Ctrl1','Ctrl2','Ctrl3','Ctrl4'), project = 'GSE174609')
pbmc

pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")
# Visualize_QC_metrics
pdf(file=paste0(outputdir,"feature_vlnplot_GSE174609.pdf"))
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
# Visualize_feature-feature_relationships
pdf(file=paste0(outputdir,"feature_scatter_GSE174609.pdf"),width=13,height=10)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()

saveRDS(pbmc, file = paste0(outputdir,"/preprocessing_GSE174609.rds"))

