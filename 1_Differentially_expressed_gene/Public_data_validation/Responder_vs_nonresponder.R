#! /usr/bin/env Rscript

# Updated by Hye-Yeon Ju on Apr 7, 2025.

##### Code Description #####
# This code performs gene expression profiling of chemotherapy responder and non-responder samples.

##### Code Example #####
# Rscript
# /home/jhy/single_cell_TNBC/code/for_github/1_Differentially_expressed_gene/Public_data_validation/Responder_vs_nonresponder.R

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)

data <- readRDS('/home/jhy/single_cell_TNBC/result/PrePostChemo_250403/UseThis_nobatch/EffTcell_dim10_res0.4.rds')
outputdir <- '/home/jhy/single_cell_TNBC/result/ResNonresChemo_250405/Post_chemo/'
outputdir <- '/home/jhy/single_cell_TNBC/result/ResNonresChemo_250405/Matched_all/'
outputdir <- '/home/jhy/single_cell_TNBC/result/ResNonresChemo_250405/Pre_chemo/'

## Only Post or Pre-treatmnet
test <- colnames(data)
include_list <- c("P022", "P011", "P020", "P008", "P013", "P025", "P018", "P023", "P024")
filtered_elements <- grep("\\.Post[^_]*_.*_b$", test, value = TRUE, perl = TRUE)
filtered_elements <- grep(paste0("^.*\\.Post[^_]*_((", paste(include_list, collapse = ")|("), ")).*_b$"), filtered_elements, value = TRUE, perl = TRUE)
data <- data[,colnames(data) %in% filtered_elements]

## Responder
test <- colnames(data)
include_list <- c("P022", "P011", "P020", "P008", "P013")
filtered_elements <- grep(paste(include_list, collapse = "|"), test, value = TRUE)
data@meta.data[filtered_elements, "orig.ident"] <- 'Responder'

## Non-responder
include_list <- c("P025", "P018", "P023", "P024")
filtered_elements <- grep(paste(include_list, collapse = "|"), test, value = TRUE)
data@meta.data[filtered_elements, "orig.ident"] <- 'Non_Responder'

## Normalization and scaling
Idents(data) <- data$orig.ident
data <- NormalizeData(data)
data <- ScaleData(data)
#library(future)
#plan(sequential)
#options(future.globals.maxSize = 32 * 1024^3)
data <- SCTransform(data)

## Extract marker list
pbmc.markers <- FindAllMarkers(object = data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
p <- subset(pbmc.markers, subset = p_val_adj <= 0.05)
write.table(data.frame(p), paste0(outputdir,'markers_log2FC0.25_adjp0.05_dim10.txt'))

pbmc.markers <- FindAllMarkers(object = data, only.pos = TRUE, min.pct = 0, logfc.threshold = 0,return.thresh = 1)
write.table(data.frame(pbmc.markers), paste0(outputdir,'markers_all.txt'))

## DotPlot
pdf(paste0(outputdir,'DotPlot.pdf'),width=8,height=5)
DotPlot(data, features=c('CST7', 'GZMA', 'ACTB', 'CCL5', 'ATP5F1E', 'IFITM1', 'ATP5MG', 'KLRB1', 'SLC9A3R1', 'RPL17', 'GZMB', 'CD74', 'HLA-B')) + theme(axis.text.x = element_text(angle=45,hjust=1)) + scale_y_discrete(limits = rev(levels(Idents(data))))
dev.off()

## VlnPlot
vp_case1 <- function(gene_signature, file_name, test_sign){
	plot_case1 <- function(signature, y_max = 5){
		VlnPlot(data, features = signature, pt.size = 0, idents = c('Responder','Non_Responder'), y.max = y_max) + stat_compare_means(comparisons = test_sign, method = "wilcox.test", label = "p.signif") + geom_boxplot(width=0.15,fill='white',color='black',outlier.shape=NA,position = position_dodge(width = 0.9),size=0.5) + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) + stat_summary(fun=mean, geom='point',size=0.5)}
        plot_list <- list()
	for (gene in gene_signature) {
		plot_list[[gene]] <- plot_case1(gene)
	}
	cowplot::plot_grid(plotlist = plot_list)
	file_name <- paste0(file_name, ".pdf")
	ggsave(file_name, width = 18, height = 18)
}
gene_sig <- c('CST7', 'GZMA', 'ACTB', 'CCL5', 'ATP5F1E', 'IFITM1', 'ATP5MG', 'KLRB1', 'SLC9A3R1', 'RPL17', 'GZMB', 'CD74', 'HLA-B')
comparisons1 <- list(c('Responder','Non_Responder'))
vp_case1(gene_signature = gene_sig, file_name = paste0(outputdir,'Vlnplot_with_boxplot.pdf'), test_sign=comparisons1)


