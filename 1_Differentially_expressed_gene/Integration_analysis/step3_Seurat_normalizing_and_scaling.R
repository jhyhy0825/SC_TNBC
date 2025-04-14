#! /usr/bin/env Rscript

# Updated by Hye-Yeon Ju on Jun 26, 2023.

##### Code Description #####
# This code covers the processes of normalizing and scaling.
# To consider cutoff, you need to check the plots from the previous results.
# nFeature_cutoff: filter in cells having less than N detected genes.
# percent.mt_cutoff: filter in cells cells in which mitochondrial protein-coding genes represented less than N% of the UMI content.

##### Code Example #####
# Rscript + resultdir + nFeature_cutoff + percent.mt_cutoff
# /home/jhy/single_cell_TNBC/code/for_github/1_Differentially_expressed_gene/Integration_analysis/step3_Seurat_normalizing_and_scaling.R /home/jhy/single_cell_TNBC/result/Seurat/ 5000 20

library(dplyr)
library(Seurat)
library(patchwork)

args = commandArgs(trailingOnly = TRUE)
outputdir <- args[1]
nFeature <- as.numeric(args[2])
per_mt <- as.numeric(args[3])

pbmc <- readRDS(file = paste0(outputdir,"/preprocessing.rds"))
DefaultAssay(pbmc) <- "RNA"
pbmc

### Quality_control ###
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < nFeature & percent.mt < per_mt)
pbmc

### Visualize_QC_metrics_after_QC ###
pdf(file=paste0(outputdir,"feature_vlnplot_after_QC.pdf"))
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
# Visualize_feature-feature_relationships_after_QC
pdf(file=paste0(outputdir,"feature_scatter_after_QC.pdf"),width=13,height=10)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()

case <- subset(pbmc, orig.ident=="8674-TK-1-Case" | orig.ident=="8674-TK-2-Case")
case[["info"]] <- "Case"
case
control <- subset(pbmc, orig.ident=="8674-TK-1-Control" | orig.ident=="8674-TK-2-Control")
control
control[["info"]] <- "Control"
pbmc <- merge(case, control)
pbmc
head(pbmc)

### Normalizing_the_data ###
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

DefaultAssay(pbmc) <- "ADT"
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
DefaultAssay(pbmc) <- "RNA"

### Identification_of highly_variable_features ###
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)
# Identify_the_10_most_highly_variable_genes
top10 <- head(VariableFeatures(pbmc), 10)
# Plot_variable_features_with_and_without_labels
pdf(file=paste0(outputdir,"variable_feature_plot.pdf"),width=13,height=10)
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
dev.off()


### Scaling_data ###
all.genes <- rownames(x = pbmc)
pbmc <- ScaleData(object = pbmc, features = all.genes)

saveRDS(pbmc, file = paste0(outputdir,"/normalizing_and_scaling.rds"))


