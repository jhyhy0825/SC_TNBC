#! /usr/bin/env Rscript

# Updated by Hye-Yeon Ju on Apr 18, 2024.

##### Code Description #####
# This code integrates internal and public datasets.
# This code performes batch correction, clustering, and visualization.

##### Code Example #####
# Rscript
# /home/jhy/single_cell_TNBC/code/integration/Harmony_with_control.R

library(dplyr)
library(Seurat)
library(patchwork)
library(harmony)

### Import RDS files of internal and public datasets ###
obj1 <- readRDS(file = '/home/jhy/single_cell_TNBC/result/Seurat/preprocessing.rds')
obj1$orig.ident <- "Vanderbilt"
obj2 <- readRDS(file = '/home/jhy/single_cell_TNBC/result/Seurat/integration/preprocessing_blood_posttreat_chemo.rds')
obj2$orig.ident <- "GSE169246_TNBC"
obj3 <- readRDS(file = '/home/jhy/single_cell_TNBC/result/Seurat/integration2/preprocessing_GSE174609.rds')
obj3$orig.ident <- "GSE174609_Control"
obj4 <- readRDS(file = '/home/jhy/single_cell_TNBC/result/Seurat/integration2/preprocessing_GSE134355.rds')
obj4$orig.ident <- "GSE134355_Control"

### Merge_datasets ###
raw_seurat_list <- c(obj1,obj2,obj3,obj4)

merged_seurat <- merge(x = raw_seurat_list[[1]],y = raw_seurat_list[2:length(raw_seurat_list)],merge.data = TRUE)
merged_seurat <- merged_seurat %>% NormalizeData() %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% ScaleData() %>% SCTransform()
#merged_seurat <- SCTransform()
merged_seurat <- RunPCA(merged_seurat, assay = "SCT")
merged_seurat

outputdir <- '/home/jhy/single_cell_TNBC/result/Seurat/integration2/result2/'

### Visualize_ElbowPlot ###
pdf(file=paste0(outputdir,"ElbowPlot.pdf"))
ElbowPlot(merged_seurat)
dev.off()

### Batch_correction_and_Clustering ###
merged_seurat
head(merged_seurat)
harmonized_seurat <- RunHarmony(merged_seurat, group.by.vars = "orig.ident", reduction = "pca", assay.use = "SCT", reduction.save = "harmony")
harmonized_seurat <- RunUMAP(harmonized_seurat, reduction = "harmony", assay = "SCT", dims = 1:10)
harmonized_seurat <- FindNeighbors(object = harmonized_seurat, reduction = "harmony", dims=1:10)
harmonized_seurat <- FindClusters(harmonized_seurat, resolution = 0.4)
harmonized_seurat

pdf(file=paste0(outputdir,"umap_dim10_res0.4.pdf"))
DimPlot(object=harmonized_seurat, reduction="umap")
dev.off()

pdf(file=paste0(outputdir,"umap_label_dim10_res0.4.pdf"))
DimPlot(object=harmonized_seurat, label=TRUE, reduction="umap", label.size=9)
dev.off()

#pdf(file=paste0(outputdir,"FeaturePlot1_dim10.pdf"), width=16, height=12)
#FeaturePlot(harmonized_seurat, features = c('Stat1','Stat3','Stat4','Stat6','Stat5a','Il6ra','Il6st','Il10ra','Il10rb','Il21r','Il27ra'))
#dev.off()

#pdf(file=paste0(outputdir,"FeaturePlot2_dim10.pdf"), width=16, height=12)
#FeaturePlot(harmonized_seurat, features = c('Pdcd1','Tcf7','Tox','Havcr2','Mki67','Sell','Rora','Il18r1','Il21r','Tnfrsf4'))
#dev.off()

saveRDS(harmonized_seurat, file = paste0(outputdir,'/Integration_dim10_res0.4.rds'))

