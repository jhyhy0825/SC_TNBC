#! /usr/bin/env Rscript

# Updated by Hye-Yeon Ju on Apr 23, 2024.

##### Code Description #####
# This code performs visualization.

##### Code Example #####
# Rscript 
# /home/jhy/single_cell_TNBC/code/for_github/1_Differentially_expressed_gene/Integration_analysis/Visualization.R

library(dplyr)
library(Seurat)
library(patchwork)
library(ggpubr)

### DotPlot ###
DotPlot(merged_seurat,feature=c('CD8A', 'CD8B', 'GZMK', 'CD4', 'FOXP3', 'CTLA4', 'TOX', 'HAVCR2', 'TNFRSF4', 'PDCD1', 'SELL', 'CCR7', 'CD27', 'CD28', 'PTPRC', 'GZMB', 'PRF1', 'CCL5',  'ITGA1', 'ITGB2', 'CD70', 'CD44', 'KLRD1', 'KLRB1', 'KLRG1', 'CX3CR1', 'CD69', 'CD58', 'CD2', 'NKG7', 'ITAM'),split.by="orig.ident", cols=c("red","blue")) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

### Heatmap ###
LogNormalize and avergeExp
new <- NormalizeData(object=data, assay='RNA')
avgexp = AverageExpression(new, assay='RNA', group.by='seurat_clusters',return.seurat=T)
test <- ScaleData(avgexp)
pdf('/home/jhy/single_cell_TNBC/result/Seurat/integration2/result2/tcell_heatmap_test4.pdf',height=8,width=10)
DoHeatmap(test,feature=c('CD8A', 'CD8B', 'GZMK', 'CD4', 'FOXP3', 'CTLA4','ICOS','TIGIT','TOX', 'HAVCR2', 'TNFRSF4', 'PDCD1', 'SELL', 'CCR7', 'CD27', 'CD28', 'PTPRC', 'GZMB', 'PRF1', 'CCL5',  'ITGA1', 'ITGB2', 'CD70', 'CD44', 'KLRD1', 'KLRB1', 'KLRG1', 'CX3CR1', 'CD69', 'CD58', 'CD2', 'NKG7', 'ITAM'))+ scale_fill_gradientn(colors = c("blue", "white", "red"))
#DoHeatmap(test,feature=c('CD8A', 'CD8B', 'GZMK', 'CD4', 'FOXP3', 'CTLA4', 'TOX', 'HAVCR2', 'TNFRSF4', 'PDCD1', 'SELL', 'CCR7', 'CD27', 'CD28', 'PTPRC', 'GZMB', 'PRF1', 'CCL5',  'ITGA1', 'ITGB2', 'CD70', 'CD44', 'KLRD1', 'KLRB1', 'KLRG1', 'CX3CR1', 'CD69', 'CD58', 'CD2', 'NKG7', 'ITAM'))+ scale_fill_gradientn(colors = c("blue", "white", "red"))
#DoHeatmap(test,feature=c('PTPRC','CD3D','CD2','NKG7','KLRF1','CD79A','MS4A1','S100A8','LYZ','EPCAM','COL1A1'))+ scale_fill_gradientn(colors = c("blue", "white", "red"))
#'NKG7','KLRF1', 'CD2', 'NCAM1', 'B3GAT1', 'CD69', 'LILRB1', 'KLRC2', 'PTPRC', 'ID2'
'CD14','S100A12','LYZ','CD163','CD68','CLEC4C','IL3RA','LILRA4','CLEC9A','CLEC10A','CD1C','FCER1A','FCGR3B','CXCR2'
## draw.lines=F
dev.off()

### BarPlot ###
cluster_ids <- Idents(merged_seurat)
orig_idents <- merged_seurat$orig.ident
grouped_data <- data.frame(cluster = factor(cluster_ids),orig_ident = orig_idents)
cluster_counts <- with(grouped_data, table(cluster, orig_ident))
print(cluster_counts)
pdf(file=paste0(outputdir,"bar_tcellall_mt10_feat4000_dim10_res0.2_casecontrol.pdf"),height=8,width=10)
library(tidyr)
cluster_data_long <- pivot_longer(cluster_data, cols = c(case, control), names_to = "group", values_to = "value")
pdf(file=paste0(outputdir,"bar_tcellall_mt10_feat4000_dim10_res0.2_casecontrol.pdf"),height=4,width=3)
ggplot(cluster_data_long, aes(x = factor(cluster), y = value, fill = group)) + geom_bar(stat = "identity", position = "stack") + labs(x = "Cluster", y = "Value") + scale_fill_manual(values = c("case" = 'coral2', "control" = 'cyan3')) + theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
markers <- FindMarkers(object = newT2, ident.1 = c(5,6,8), ident.2=c(3,7,10), min.pct = 0.25, logfc.threshold = 0.25)
cluster_data <- as.data.frame.matrix(cluster_counts)
cluster_data$cluster <- rownames(cluster_data)
colnames(cluster_data) <- c("case", "control", "cluster")
cluster_data$total <- rowSums(cluster_data[, c("case", "control")])
cluster_data$case_prop <- cluster_data$case / cluster_data$total
cluster_data$control_prop <- cluster_data$control / cluster_data$total
cluster_data_long <- pivot_longer(cluster_data, cols = c(case_prop, control_prop),names_to = "group", values_to = "proportion")
print(head(cluster_data_long))
cluster_data_long$cluster <- factor(cluster_data_long$cluster, levels = as.character(sort(as.numeric(unique(cluster_data_long$cluster)))))
pdf('/home/jhy/single_cell_TNBC/result/bar_240823/myeloid_barplot_casecontrol.pdf', height=4, width=3)
ggplot(cluster_data_long, aes(x = cluster, y = proportion, fill = group)) +geom_bar(stat = "identity", position = "stack") +labs(x = "Cluster", y = "Proportion", fill = "Group") +scale_fill_manual(values = c("case_prop" = 'coral2', "control_prop" = 'cyan3')) +theme_minimal() +theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

### VlnPlot ###
genes <- c('PTPRCAP','HOPX','GZMK','CTSW','KLRG1','PRF1','NKG7','JUNB','CCL4','GNLY')
genes <- c('HOPX','CST7','KLF6','GZMA','ACTB','CCL5','ATP5F1E','IFITM1','ATP5MG','KLRB1','SLC9A3R1','RPL17','GZMB','CD74','HLA-B')
vp_case1 <- function(gene_signature, file_name,test_sign){
	        plot_case1 <- function(signature, y_max = 5){
			                VlnPlot(pbmc, features = signature, pt.size = 0, group.by = "orig.ident", y.max = y_max) + stat_compare_means(comparisons=test_sign, method = "wilcox.test", label = "p.signif")}
plot_list <- list()
for (gene in gene_signature) {
	        plot_list[[gene]] <- plot_case1(gene)
}
cowplot::plot_grid(plotlist = plot_list)
file_name <- paste0(file_name, ".pdf")
ggsave(file_name, width = 10, height = 9)
}
gene_sig <- genes
comparisons1 <- list(c("case","control"))
vp_case1(gene_signature = gene_sig, file_name = '/home/jhy/single_cell_TNBC/result/Seurat/integration2/result2/240705_vlnplot_top10genes.pdf', test_sign = comparisons1)

