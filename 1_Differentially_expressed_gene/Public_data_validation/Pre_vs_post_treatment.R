#! /usr/bin/env Rscript

# Updated by Hye-Yeon Ju on Apr 7, 2025.

##### Code Description #####
# This code performs gene expression profiling of pre- and post-chemotherapy treated samples.

##### Code Example #####
# Rscript
# /home/jhy/single_cell_TNBC/code/for_github/1_Differentially_expressed_gene/Public_data_validation/Pre_vs_post_treatment.R 

library(Seurat)
library(dplyr)
library(harmony)
library(ggpubr)

## Preprocessing이 끝난 GSE169246 데이터 (전체)
data <- readRDS('/data1/analysis/BC_INU_VU_KU/03_TNBC_metastasis_scRNAseq/BCG_analysis/05_results/Seurat/integration/preprocessing.rds')

## Filter in Post samples
test <- colnames(data)
include_list <- c("P022", "P011", "P020", "P008", "P013", "P025", "P018", "P023", "P024", "P003", "P028")
filtered_elements <- grep("\\.Post[^_]*_.*_b$", test, value = TRUE, perl = TRUE)
filtered_elements <- grep(paste0("^.*\\.Post[^_]*_((", paste(include_list, collapse = ")|("), ")).*_b$"), filtered_elements, value = TRUE, perl = TRUE)
new <- data[,colnames(data) %in% filtered_elements]
test2 <- colnames(new)
extracted_strings <- sapply(strsplit(test2, "_"), function(x) x[2])
saveRDS(new,'/home/jhy/single_cell_TNBC/result/PrePostChemo_250403/post_chemo_blood.rds')

## Filter in Pre samples
#For All
#test <- colnames(data)
#include_list <- c("P022", "P011", "P020", "P008", "P013", "P025", "P018", "P023", "P024", "P003", "P028")
#For Matched
test <- colnames(data)
include_list <- c("P022", "P011", "P020", "P008", "P013", "P025", "P018", "P023", "P024", "P003")
filtered_elements <- grep("\\.Pre[^_]*_.*_b$", test, value = TRUE, perl = TRUE)
filtered_elements <- grep(paste0("^.*\\.Pre[^_]*_((", paste(include_list, collapse = ")|("), ")).*_b$"), filtered_elements, value = TRUE, perl = TRUE)
new <- data[,colnames(data) %in% filtered_elements]
test2 <- colnames(new)
extracted_strings <- sapply(strsplit(test2, "_"), function(x) x[2])
saveRDS(new,'/home/jhy/single_cell_TNBC/result/PrePostChemo_250403/pre_chemo_blood.rds')

## Integration
post <- readRDS('/home/jhy/single_cell_TNBC/result/PrePostChemo_250403/post_chemo_blood.rds')
post$orig.ident <- 'Post_chemo'
pre <- readRDS('/home/jhy/single_cell_TNBC/result/PrePostChemo_250403/pre_chemo_blood.rds')
pre$orig.ident <- 'Pre_chemo'

raw_seurat_list <- c(post,pre)
merged_seurat <- merge(x = raw_seurat_list[[1]],y = raw_seurat_list[2:length(raw_seurat_list)],merge.data = TRUE)
merged_seurat <- merged_seurat %>% NormalizeData() %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% ScaleData() %>% SCTransform()
##만약 여기서 SCTransform 에러가 난다면##
#library(future)
#plan(sequential)
#options(future.globals.maxSize = 32 * 1024^3)
#merged_seurat <- SCTransform(merged_seurat)
##########################################
merged_seurat <- RunPCA(merged_seurat, assay = "SCT")
merged_seurat
outputdir <- '/home/jhy/single_cell_TNBC/result/PrePostChemo_250403/'
pdf(file=paste0(outputdir,"ElbowPlot.pdf"))
ElbowPlot(merged_seurat)
dev.off()
#saveRDS(merged_seurat,'/home/jhy/single_cell_TNBC/result/PrePostChemo_250403/merged_data.rds')
##Dicide the number of PCs uisng ElbowPlot result

## Batch correction
#harmonized_seurat <- RunHarmony(merged_seurat, group.by.vars = "orig.ident", reduction = "pca", assay.use = "SCT", reduction.save = "harmony")
#harmonized_seurat <- RunUMAP(harmonized_seurat, reduction = "harmony", assay = "SCT", dims = 1:10)
#harmonized_seurat <- FindNeighbors(object = harmonized_seurat, reduction = "harmony", dims=1:10)
#harmonized_seurat <- FindClusters(harmonized_seurat, resolution = 0.4)
#harmonized_seurat

outputdir <- '/home/jhy/single_cell_TNBC/result/PrePostChemo_250403/UseThis_nobatch/'

## No Batch corection (Ignore name)
harmonized_seurat <- FindNeighbors(object = merged_seurat, dims = 1:10)
harmonized_seurat <- FindClusters(object = harmonized_seurat, resolution = 0.4)
harmonized_seurat <- RunUMAP(object = harmonized_seurat, dims = 1:10)

pdf(file=paste0(outputdir,"umap_dim10_res0.4.pdf"))
DimPlot(object=harmonized_seurat, reduction="umap")
dev.off()

pdf(file=paste0(outputdir,"umap_label_dim10_res0.4.pdf"))
DimPlot(object=harmonized_seurat, label=TRUE, reduction="umap", label.size=9)
dev.off()

saveRDS(harmonized_seurat, file = paste0(outputdir,'/pre_and_post_dim10_res0.4.rds'))


## Cell type annotation
library(scCATCH)
Seurat_obj <- harmonized_seurat

# Main Cell type
gene = c('PTPRC', 'CD2', 'CD3D', 'CD3E', 'CD5', 'CD6', 'CD7', 'CD4', 'CD8A', 'CD8B', 'CD27', 'CD28', 'CD3G', 'UCHL1', 'CCR7', 'IL7R', 'IL2RA', 'BTLA', 'SELL', 'CTLA4', 'PDCD1', 'HAVCR2', 'CD247', 'CD14', 'CD163', 'CD68', 'CSF1R', 'FCGR3A', 'ITGAM', 'LYZ', 'MRC1', 'ITGAX', 'CD83', 'ITGA2B', 'GP1BA', 'ITGAM', 'ITGAX', 'IL3RA', 'FUT4', 'FCGR3A', 'CD33', 'PTPRC', 'CEACAM8', 'CASK', 'CD14', 'BLNK', 'CD19', 'CD79A', 'CD79B', 'MS4A1', 'CD37', 'MME', 'TNFSF13B', 'CD28', 'CD38', 'TFRC', 'CD86', 'CR2', 'CD24', 'CD27', 'PTPRC', 'CD5', 'IGHD', 'IGHM', 'ZAP70', 'CD74' , 'CCL3', 'CD160', 'CD247', 'GNLY', 'GZMB', 'NKG7', 'FCGR3B', 'KLRB1', 'KLRC1', 'KLRD1', 'KLRF1', 'KLRK1', 'NCAM1', 'ITGAM', 'FCGR3A', 'CD19', 'CD2', 'CD3D', 'CD3E', 'PTPRC', 'ITGA1', 'CXCR6', 'CLEC4C', 'IL3RA', 'NRP1', 'CLEC9A', 'ITGAX', 'CD80', 'CD83', 'CD86', 'CD34', 'KRT19', 'CD24', 'CDH1', 'CLDN1', 'EPCAM', 'KRT14', 'KRT7', 'MUC1', 'BMX', 'CD93', 'ECE1', 'SELE')
celltype = c('T cell', 'T cell', 'T cell', 'T cell', 'T cell', 'T cell', 'T cell', 'T cell', 'T cell', 'T cell', 'T cell', 'T cell', 'T cell', 'T cell', 'T cell', 'T cell', 'T cell', 'T cell', 'T cell', 'T cell', 'T cell', 'T cell', 'T cell', 'Macrophage', 'Macrophage', 'Macrophage', 'Macrophage', 'Macrophage', 'Macrophage', 'Macrophage', 'Macrophage', 'Macrophage', 'Macrophage', 'Myeloid cell', 'Myeloid cell', 'Myeloid cell', 'Myeloid cell', 'Myeloid cell', 'Myeloid cell', 'Myeloid cell', 'Myeloid cell', 'Myeloid cell', 'Myeloid cell', 'Myeloid cell', 'Myeloid cell', 'B cell', 'B cell', 'B cell', 'B cell', 'B cell', 'B cell', 'B cell', 'B cell', 'B cell', 'B cell', 'B cell', 'B cell', 'B cell', 'B cell', 'B cell', 'B cell', 'B cell', 'B cell', 'B cell', 'B cell', 'B cell', 'Natural killer cell', 'Natural killer cell', 'Natural killer cell', 'Natural killer cell', 'Natural killer cell', 'Natural killer cell', 'Natural killer cell', 'Natural killer cell', 'Natural killer cell', 'Natural killer cell', 'Natural killer cell', 'Natural killer cell', 'Natural killer cell', 'Natural killer cell', 'Natural killer cell', 'Natural killer cell', 'Natural killer cell', 'Natural killer cell', 'Natural killer cell', 'Natural killer cell', 'Natural killer cell', 'Natural killer cell', 'dendritic cell', 'dendritic cell', 'dendritic cell', 'dendritic cell', 'dendritic cell', 'dendritic cell', 'dendritic cell', 'dendritic cell', 'dendritic cell', 'Epithelial cell', 'Epithelial cell', 'Epithelial cell', 'Epithelial cell', 'Epithelial cell', 'Epithelial cell', 'Epithelial cell', 'Epithelial cell', 'Endothelial cell', 'Endothelial cell', 'Endothelial cell', 'Endothelial cell')


custom_marker <- data.frame(species="Human", tissue="Breast", cancer="Breast Cancer", condition=NA,subtype1=NA,subtype2=NA,subtype3=NA,celltype=celltype,gene=gene,resource=NA,pmid=NA,stringsAsFactors=FALSE)

obj <- createscCATCH(data = Seurat_obj[["SCT"]]@data, cluster = as.character(Idents(Seurat_obj)))
obj <- findmarkergene(object=obj, if_use_custom_marker = TRUE, species="Human", marker=custom_marker, tissue="Breast", use_method="2")
obj <- findcelltype(obj)
obj@celltype

pdf(paste0(outputdir,'FeaturePlot_main_cell.pdf'),width=12,height=8)
FeaturePlot(Seurat_obj, features=c('PTPRC','CD3D','CD2','NKG7','KLRF1','CD79A','MS4A1','S100A8','LYZ','EPCAM','COL1A1'))
dev.off()

pdf(paste0(outputdir,'DotPlot_main_cell.pdf'),width=12,height=8)
DotPlot(Seurat_obj, features=c('PTPRC','CD3D','CD2','NKG7','KLRF1','CD79A','MS4A1','S100A8','LYZ','EPCAM','COL1A1'))
dev.off()

## Extract marker list
pbmc.markers <- FindAllMarkers(object = Seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
p <- subset(pbmc.markers, subset = p_val_adj <= 0.05)
write.table(data.frame(p), paste0(outputdir,'markers_log2FC0.25_adjp0.05_dim10.txt'))

## Sub-clustering
Tcell <- subset(Seurat_obj, seurat_clusters %in% c(0,1,3,8,9,13,16))
pbmc <- RunPCA(object = Tcell, features = VariableFeatures(object = Tcell))
pdf(file=paste0(outputdir,"Tcell_ElbowPlot.pdf"))
ElbowPlot(pbmc)
dev.off()

pbmc <- FindNeighbors(object = pbmc, dims = 1:10)
pbmc <- FindClusters(object = pbmc, resolution = 0.4)
pbmc <- RunUMAP(object = pbmc, dims = 1:10)
pdf(file=paste0(outputdir,"Tcell_umap.pdf"))
DimPlot(object = pbmc, reduction = "umap")
dev.off()
pdf(file=paste0(outputdir,"Tcell_umap_label.pdf"))
DimPlot(object = pbmc, reduction = "umap", label=TRUE, label.size=9)
dev.off()
pdf(file=paste0(outputdir,"Tcell_umap_label_pt.pdf"))
DimPlot(object = pbmc, reduction = "umap", label=TRUE, label.size=9, pt.size=0.2)
dev.off()
pdf(file=paste0(outputdir,"Tcell_umap_patient.pdf"))
DimPlot(object = pbmc, reduction = "umap",split.by='orig.ident')
dev.off()
Idents(pbmc) <- pbmc$orig.ident
pdf(file=paste0(outputdir,"Tcell_umap_patient2.pdf"))
DimPlot(object = pbmc, reduction = "umap")
dev.off()
Idents(pbmc) <- pbmc$seurat_clusters
pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
p <- subset(pbmc.markers, subset = p_val_adj <= 0.05)
write.table(data.frame(p), paste0(outputdir,'Tcell_markers_log2FC0.25_adjp0.05.txt'))
saveRDS(pbmc, paste0(outputdir,'Tcell_pre_and_post_dim10_res0.4.rds'))

## Cell type annotation
gene = c('CD8A', 'CD8B', 'GZMK', 'TNFRSF8', 'NCAM1', 'NKG7', 'ITGAM', 'CD3', 'TOX', 'HAVCR2', 'TNFRSF4', 'PDCD1' , 'CD4', 'CTLA4', 'FOXP3', 'IL2RA', 'CCR6', 'LTB', 'SLC2A1', 'CD3', 'CCR10', 'CD52', 'CMTM7', 'FOXP3', 'CTLA4', 'FOXP3', 'IL2RA', 'ENTPD1', 'CD4', 'IL7R', 'ITGA4', 'VPS53', 'TNFRSF18', 'CD101', 'IL2RB', 'FCRL3', 'TIGIT', 'IKZF2', 'CD3', 'IL7R', 'CCR7', 'SELL', 'CD27', 'CD28' ,'CCR7', 'PTPRC', 'SELL', 'CD44', 'KLRG1' ,'GZMB', 'PRF1', 'CD8A', 'CD8B', 'CCL5', 'ITGA1', 'ITGB2', 'CTLA4', 'CD28', 'CD27', 'CD70', 'CD44', 'KLRD1', 'KLRB1', 'KLRG1', 'CX3CR1', 'CCR7', 'CD69', 'CD58', 'CD2', 'NKG7', 'ITAM')
celltype = c('cytotoxic T cell', 'cytotoxic T cell', 'cytotoxic T cell', 'cytotoxic T cell', 'cytotoxic T cell', 'cytotoxic T cell', 'cytotoxic T cell', 'cytotoxic T cell' , 'Exhausted T cell', 'Exhausted T cell', 'Exhausted T cell', 'Exhausted T cell', 'CD4 T cell', 'CD4 T cell', 'CD4 T cell', 'CD4 T cell', 'CD4 T cell', 'CD4 T cell', 'CD4 T cell', 'CD4 T cell', 'regulatory T cell', 'regulatory T cell', 'regulatory T cell','regulatory T cell', 'regulatory T cell', 'regulatory T cell', 'regulatory T cell', 'regulatory T cell', 'regulatory T cell', 'regulatory T cell', 'regulatory T cell', 'regulatory T cell', 'regulatory T cell', 'regulatory T cell', 'regulatory T cell', 'regulatory T cell', 'regulatory T cell', 'regulatory T cell', 'regulatory T cell', 'Central Memory T cell', 'Central Memory T cell', 'Central Memory T cell', 'Central Memory T cell', 'Central Memory T cell', 'Effector Memory T cell', 'Effector Memory T cell', 'Effector Memory T cell', 'Effector Memory T cell', 'Effector Memory T cell', 'Terminally Differentiated T cell', 'Terminally Differentiated T cell', 'Terminally Differentiated T cell', 'Terminally Differentiated T cell', 'Terminally Differentiated T cell', 'Terminally Differentiated T cell', 'Terminally Differentiated T cell', 'Terminally Differentiated T cell', 'Terminally Differentiated T cell', 'Terminally Differentiated T cell', 'Terminally Differentiated T cell', 'Terminally Differentiated T cell', 'Terminally Differentiated T cell', 'Terminally Differentiated T cell', 'Terminally Differentiated T cell', 'Terminally Differentiated T cell', 'Terminally Differentiated T cell', 'Terminally Differentiated T cell', 'Terminally Differentiated T cell', 'Terminally Differentiated T cell', 'Terminally Differentiated T cell', 'Terminally Differentiated T cell')
custom_marker <- data.frame(species="Human", tissue="Breast", cancer="Breast Cancer", condition=NA,subtype1=NA,subtype2=NA,subtype3=NA,celltype=celltype,gene=gene,resource=NA,pmid=NA,stringsAsFactors=FALSE)
obj <- createscCATCH(data = pbmc[["SCT"]]@data, cluster = as.character(Idents(pbmc)))
obj <- findmarkergene(object=obj, if_use_custom_marker = TRUE, species="Human", marker=custom_marker, tissue="Breast", use_method="2")
obj <- findcelltype(obj)
obj@celltype

library(ggplot2)
pdf(paste0(outputdir,'DotPlot_Tcell.pdf'),width=10,height=4)
DotPlot(pbmc, features = c('CD8A','CD8B','GZMK','CD4','FOXP3','CTLA4','TOX','HAVCR2','TNFRSF4','PDCD1','SELL','CCR7','CD27','CD28','PTPRC','GZMB','PRF1','CCL5','ITGA1','ITGB2','CD70','CD44','KLRD1','KLRB1','KLRG1','CX3CR1','CD69','CD58','CD2','NKG7')) + theme(axis.text.x = element_text(angle=45,hjust=1))
dev.off()

## Effector T cell
Eff <- subset(pbmc,seurat_clusters %in% c(1,4,6,10,11,14,16,17,18))
Idents(Eff) <- Eff$orig.idens

pbmc.markers <- FindAllMarkers(object = Eff, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
p <- subset(pbmc.markers, subset = p_val_adj <= 0.05)
write.table(data.frame(p), paste0(outputdir,'EffTcell_markers_log2FC0.25_adjp0.05.txt'))

## DotPlot
pdf(paste0(outputdir,'DotPlot_EffTcell.pdf'),width=8,height=5)
DotPlot(Eff, features = c('CST7', 'GZMA', 'ACTB', 'CCL5', 'ATP5F1E', 'IFITM1', 'ATP5MG', 'KLRB1', 'SLC9A3R1', 'RPL17', 'GZMB', 'CD74', 'HLA-B')) + theme(axis.text.x = element_text(angle=45,hjust=1))
dev.off()

## VlnPlot
pdf(paste0(outputdir,'VlnPlot_EffTcell_1.pdf'),width=10,height=10)
#VlnPlot(Eff, features = c('CST7', 'GZMA', 'ACTB', 'CCL5', 'ATP5F1E', 'IFITM1', 'ATP5MG', 'KLRB1', 'SLC9A3R1', 'RPL17', 'GZMB', 'CD74', 'HLA-B'))
VlnPlot(Eff, features = c('CST7', 'GZMA', 'ACTB', 'CCL5', 'ATP5F1E', 'IFITM1'),pt.size=0)
VlnPlot(Eff, features = c('ATP5MG', 'KLRB1', 'SLC9A3R1', 'RPL17', 'GZMB', 'CD74', 'HLA-B'),pt.size=0)
dev.off()

## Heatmap
pdf(paste0(outputdir,'Heatmap_EffTcell.pdf'),width=10,height=10)
DoHeatmap(Eff, features = c('CST7', 'GZMA', 'ACTB', 'CCL5', 'ATP5F1E', 'IFITM1', 'ATP5MG', 'KLRB1', 'SLC9A3R1', 'RPL17', 'GZMB', 'CD74', 'HLA-B'))
dev.off()

## VlnPlot with boxplot
vp_case1 <- function(gene_signature, file_name, test_sign){
	plot_case1 <- function(signature, y_max = 5){
		VlnPlot(Eff, features = signature, pt.size = 0, idents = c('Post_chemo','Pre_chemo'), y.max = y_max) + stat_compare_means(comparisons = test_sign, method = "wilcox.test", label = "p.signif") + geom_boxplot(width=0.15,fill='white',color='black',outlier.shape=NA,position = position_dodge(width = 0.9),size=0.5) + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) + stat_summary(fun=mean, geom='point',size=0.5)}
	plot_list <- list()
	for (gene in gene_signature) {
		plot_list[[gene]] <- plot_case1(gene)
	}
	cowplot::plot_grid(plotlist = plot_list)
	file_name <- paste0(file_name, ".pdf")
	ggsave(file_name, width = 18, height = 18)
}

vp_case1 <- function(gene_signature, file_name, test_sign){
	plot_case1 <- function(signature, y_max = 5){
		df <- FetchData(Eff, vars = c(signature, "ident"))
		colnames(df) <- c("expression", "group")
		ggplot(df, aes(x = group, y = expression, fill = group)) + geom_boxplot(width = 0.3, outlier.shape = NA, color = "black", fill = "white", size = 0.5) + stat_summary(fun = mean, geom = "point", shape = 21, size = 1.2, fill = "black") + stat_compare_means(comparisons = test_sign, method = "wilcox.test", label = "p.signif") + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)), limits = c(0, y_max)) + theme_classic() + labs(title = signature, y = "Expression", x = NULL) + theme(axis.text.x = element_text(size = 10), plot.title = element_text(hjust = 0.5))
	}
	plot_list <- list()
	for (gene in gene_signature) {
		plot_list[[gene]] <- plot_case1(gene)
	}
	cowplot::plot_grid(plotlist = plot_list)
	file_name <- paste0(file_name, ".pdf")
	ggsave(file_name, width = 18, height = 18)
}

gene_sig <- c('CST7', 'GZMA', 'ACTB', 'CCL5', 'ATP5F1E', 'IFITM1', 'ATP5MG', 'KLRB1', 'SLC9A3R1', 'RPL17', 'GZMB', 'CD74', 'HLA-B')
comparisons1 <- list(c('Post_chemo','Pre_chemo'))
vp_case1(gene_signature = gene_sig, file_name = paste0(outputdir,'vln_withbox_test_2.pdf'), test_sign=comparisons1)


