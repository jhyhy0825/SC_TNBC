#! /usr/bin/env Rscript

# Updated by Hye-Yeon Ju on Jun 26, 2023.

##### Code Description #####
# This code covers the processes of preprocessing (merge, quality control, and visualization).

##### Code Example #####
# Rscript + cellranger_outputdir(Case1) + cellranger_outputdir(Case2) + cellranger_outputdir(Control1) + cellranger_outputdir(Control2)  + resultdir 
# /home/jhy/single_cell_TNBC/code/for_github/1_Differentially_expressed_gene/Integration_analysis/step2_Seurat_preprocessing.R /home/jhy/single_cell_TNBC/result/GEX/8674-TK-1-Case/outs/filtered_feature_bc_matrix/ /home/jhy/single_cell_TNBC/result/GEX/8674-TK-2-Case/outs/filtered_feature_bc_matrix/ /home/jhy/single_cell_TNBC/result/GEX/8674-TK-1-Control/outs/filtered_feature_bc_matrix/ /home/jhy/single_cell_TNBC/result/GEX/8674-TK-2-Control/outs/filtered_feature_bc_matrix/ /home/jhy/single_cell_TNBC/result/Seurat/

library(dplyr)
library(Seurat)
library(patchwork)

args = commandArgs(trailingOnly = TRUE)
inputdir_Case1 <- args[1]
inputdir_Case2 <- args[2]
inputdir_Control1 <- args[3]
inputdir_Control2 <- args[4]
outputdir <- args[5]

## Import_input_data
# Case1
Case1data <- Read10X(data.dir = paste0(inputdir_Case1))
Case1 <- CreateSeuratObject(counts = Case1data$'Gene Expression', project = "8674-TK-1-Case", assay='RNA', min.cells = 3, min.features = 200)
Case1
Case1[["ADT"]] <- CreateAssayObject(counts=Case1data$'Antibody Capture')
Case1
# Case2
Case2data <- Read10X(data.dir = paste0(inputdir_Case2))
Case2 <- CreateSeuratObject(counts = Case2data$'Gene Expression', project = "8674-TK-2-Case", assay='RNA', min.cells = 3, min.features = 200)
Case2
Case2[["ADT"]] <- CreateAssayObject(counts=Case2data$'Antibody Capture')
Case2
# Control1
Control1data <- Read10X(data.dir = paste0(inputdir_Control1))
Control1 <- CreateSeuratObject(counts = Control1data$'Gene Expression', project = "8674-TK-1-Control", assay='RNA', min.cells = 3, min.features = 200)
Control1
Control1[["ADT"]] <- CreateAssayObject(counts=Control1data$'Antibody Capture')
Control1
# Control2
Control2data <- Read10X(data.dir = paste0(inputdir_Control2))
Control2 <- CreateSeuratObject(counts = Control2data$'Gene Expression', project = "8674-TK-2-Control", assay='RNA', min.cells = 3, min.features = 200)
Control2
Control2[["ADT"]] <- CreateAssayObject(counts=Control2data$'Antibody Capture')
Control2

## Merge_datasets
Result <- merge(Case1, y=c(Case2,Control1,Control2), add.cell.ids=c("8674-TK-1-Case","8674-TK-2-Case","8674-TK-1-Control","8674-TK-2-Control"), project="8674-TK")
Result
head(Result)
table(Result$orig.ident)

data <- Result
data[["percent.mt"]] <- PercentageFeatureSet(object = data, pattern = "^MT-")
data

## Visualization
pdf(file=paste0(outputdir,"feature_vlnplot.pdf"))
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

pdf(file=paste0(outputdir,"feature_scatter.pdf"),width=13,height=10)
plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()

## Save_RDS_file
saveRDS(data, file = paste0(outputdir,"/preprocessing.rds"))


