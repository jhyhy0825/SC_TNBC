#! /usr/bin/env Rscript

# Updated by Hye-Yeon Ju on Nov 23, 2023.

##### Code Description #####
# This code imports public data (GSE169246, mTNBC).

##### Code Example #####
# Rscript
# /home/jhy/single_cell_TNBC/code/for_github/1_Differentially_expressed_gene/Integration_analysis/integration/For_GSE169246.R

library(dplyr)
library(Seurat)
library(patchwork)

outputdir <- '/home/jhy/single_cell_TNBC/result/Seurat/integration/'

data <- Read10X(data.dir = '/home/jhy/single_cell_TNBC/files/GSE169246/',gene.column=1)
pbmc <- CreateSeuratObject(counts = data, project = 'GSE169246', min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")
# Visualize_QC_metrics
pdf(file=paste0(outputdir,"feature_vlnplot.pdf"))
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
# Visualize_feature-feature_relationships
pdf(file=paste0(outputdir,"feature_scatter.pdf"),width=13,height=10)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()

saveRDS(pbmc, file = paste0(outputdir,"/preprocessing.rds"))


 
