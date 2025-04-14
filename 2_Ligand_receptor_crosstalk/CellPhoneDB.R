#! /usr/bin/env Rscript

# Updated by Hye-Yeon Ju on Oct 14, 2024.

##### Code Description #####
# This code performs ligand-receptor crosstalk analysis using CellPhoneDB.

##### Code Example #####
# Rscript
# /home/jhy/single_cell_TNBC/code/for_github/2_Ligand_receptor_crosstalk/CellPhoneDB.R

library(Seurat)

### Import internal data
data <- readRDS('/data1/analysis/BC_INU_VU_KU/03_TNBC_metastasis_scRNAseq/BCG_analysis/05_results/Seurat/integration2/result2/Integration_mt10_feat4000_dim10_res0.4_with_case_control.rds')

### Specify clusters based on the cell type annotation result
data@meta.data$cell_type[data@meta.data$seurat_clusters %in% c(0,1,3,6,10,11,14,16)] <- "T cell"
data@meta.data$cell_type[data@meta.data$seurat_clusters %in% c(5,12,15)] <- "B cell"
data@meta.data$cell_type[data@meta.data$seurat_clusters %in% c(2)] <- "NK cell"
data@meta.data$cell_type[data@meta.data$seurat_clusters %in% c(4,7,8,9,13,17)] <- "Myeloid cell"
case <- subset(data, orig.ident == 'case')
control <- subset(data, orig.ident == 'control')

data <- case
meta <- data.frame(data$cell_type)
meta$Cells <- rownames(meta)
colnames(meta) <- c('cell_type','Cells')
rownames(meta) <- NULL
meta <- meta %>% select(Cells,cell_type)
#meta$Cells <- gsub("-", ".", meta$Cells)
head(meta)
dim(meta)
write.table(meta,'/home/jhy/single_cell_TNBC/result/LR_241004/input/case_meta.txt',sep='\t',quote=F,row.names=F)

count_raw <- data@assays$SCT@counts[,colnames(data)]
head(count_raw)
dim(count_raw)
#write.table(count_raw,'/home/jhy/single_cell_TNBC/result/LR_241004/input/case_count.txt',sep='\t',quote=F)
#library(SeuratDisk)
#library(zellkonverter)
#library(SingleCellExperiment)
#sce <- SingleCellExperiment(assays = list(counts = as.matrix(count_raw)))
#zellkonverter::writeH5AD(sce,'/home/jhy/single_cell_TNBC/result/LR_241004/input/case_count.h5ad')
library(data.table)
count_raw <- as.data.frame(count_raw)
fwrite(count_raw, file = '/home/jhy/single_cell_TNBC/result/LR_241004/input/case_count.txt',sep='\t',quote=F,row.names=TRUE,col.names=TRUE)


### activate anaconda env cpdb ###
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method
cpdb_results = cpdb_statistical_analysis_method.call(cpdb_file_path='/home/jhy/single_cell_TNBC/code/integration/LR/cellphonedb.zip',meta_file_path='/home/jhy/single_cell_TNBC/result/LR_241004/input/case_meta.txt',counts_file_path='/home/jhy/single_cell_TNBC/result/LR_241004/input/case_count.txt',counts_data = 'hgnc_symbol',output_path = '/home/jhy/single_cell_TNBC/result/LR_241004/output/case/')
#from cellphonedb.src.core.methods import cpdb_degs_analysis_method
#cpdb_results = cpdb_degs_analysis_method.call(cpdb_file_path = '/home/jhy/single_cell_TNBC/code/integration/LR/cellphonedb.zip', meta_file_path = '/home/jhy/single_cell_TNBC/result/LR_241004/input/method2_allcell_thr0.2_p0.05/case_meta.txt', counts_file_path = '/home/jhy/single_cell_TNBC/result/LR_241004/input/method2_allcell_thr0.2_p0.05/case_count.txt', degs_file_path = '/home/jhy/single_cell_TNBC/result/LR_241004/input/method3_allcell_thr0.2_p0.05/Allcell_subset_markers_log2FC0.25_padj0.05.txt', counts_data = 'hgnc_symbol', output_path = '/home/jhy/single_cell_TNBC/result/LR_241004/output/method3_allcell_thr0.2_p0.05/case/')
cpdb_results = cpdb_degs_analysis_method.call(cpdb_file_path = '/home/jhy/single_cell_TNBC/code/integration/LR/cellphonedb.zip', meta_file_path = '/home/jhy/single_cell_TNBC/result/LR_241004/input/method2_allcell_thr0.2_p0.05/control_meta.txt', counts_file_path = '/home/jhy/single_cell_TNBC/result/LR_241004/input/method2_allcell_thr0.2_p0.05/control_count.txt', degs_file_path = '/home/jhy/single_cell_TNBC/result/LR_241004/input/method3_allcell_thr0.2_p0.05/Allcell_subset_markers_log2FC0.25_padj0.05.txt', counts_data = 'hgnc_symbol', output_path = '/home/jhy/single_cell_TNBC/result/LR_241004/output/method3_allcell_thr0.2_p0.05/control/')
cpdb_results = cpdb_degs_analysis_method.call(cpdb_file_path = '/home/jhy/single_cell_TNBC/code/integration/LR/cellphonedb.zip', meta_file_path = '/home/jhy/single_cell_TNBC/result/LR_241004/input/method2_tcell_thr0.2_p0.05/case_meta.txt', counts_file_path = '/home/jhy/single_cell_TNBC/result/LR_241004/input/method2_tcell_thr0.2_p0.05/case_count.txt', degs_file_path = '/home/jhy/single_cell_TNBC/result/LR_241004/input/method3_tcell_thr0.2_p0.05/Tcell_subset_markers_log2FC0.25_padj0.05.txt', counts_data = 'hgnc_symbol', output_path = '/home/jhy/single_cell_TNBC/result/LR_241004/output/method3_tcell_thr0.2_p0.05/case/')
cpdb_results = cpdb_degs_analysis_method.call(cpdb_file_path = '/home/jhy/single_cell_TNBC/code/integration/LR/cellphonedb.zip', meta_file_path = '/home/jhy/single_cell_TNBC/result/LR_241004/input/method2_tcell_thr0.2_p0.05/control_meta.txt', counts_file_path = '/home/jhy/single_cell_TNBC/result/LR_241004/input/method2_tcell_thr0.2_p0.05/control_count.txt', degs_file_path = '/home/jhy/single_cell_TNBC/result/LR_241004/input/method3_tcell_thr0.2_p0.05/Tcell_subset_markers_log2FC0.25_padj0.05.txt', counts_data = 'hgnc_symbol', output_path = '/home/jhy/single_cell_TNBC/result/LR_241004/output/method3_tcell_thr0.2_p0.05/control/')

import pandas as pd
df = pd.read_csv('/home/jhy/single_cell_TNBC/result/LR_241004/output/case/statistical_analysis_pvalues_10_07_2024_110111.txt',delimiter='\t')
df.head()
numeric_columns = df.columns[14:]
df[numeric_columns] = df[numeric_columns].apply(pd.to_numeric, errors='coerce')
filtered_df = df[(df[numeric_columns] <= 0.05).any(axis=1)]
ligand_receptor_df = filtered_df[filtered_df['directionality'] == 'Ligand-Receptor']
ligand_receptor_df.head()
ligand_receptor_df.to_csv('/home/jhy/single_cell_TNBC/result/LR_241004/output/case/case_p0.05_LR.txt',sep='\t',index=False)

###### p value 파일을  필터링 한 다음에 결과에서 ligand-receptor랑 simple만 선택해서 필터링한 파일을 사용 #####

data <- read.table('/home/jhy/single_cell_TNBC/result/LR_241004/output/case/case_p0.05_LR_simple.txt',sep='\t',header=T)
#cells <- colnames(data)[14:29]
cells <- colnames(data)[14:22]
cells <- lapply(cells, function(x) gsub('cell.','cell - ',x))
cells <- lapply(cells, function(x) gsub('\\.', ' ', x))
colnames(data) <- c(colnames(data)[1:13],cells)
create_combinations <- function(list_a) {
        result <- character()
        for (i in seq_along(list_a)) {
                for (j in seq_along(list_a)) {
                        result <- c(result, paste0(list_a[i], " - ", list_a[j]))
                }
        }
        return(result)
}
#a <- c('B cell','Myeloid cell','NK cell','T cell')
a <- c('Central memory T cell','Effector T cell','Regulatory T cell')
b <- create_combinations(a)
colname <- c(colnames(data)[1:13],b)
new <- data[,colname]
library(reshape2)
#test <- t(new[,14:29])
test <- t(new[,14:22])
colnames(test) <- new$interacting_pair
melted_df <- melt(test, id.vars = rownames(test))
colnames(melted_df) <- c("cell_cell", "gene_gene", "pval")
head(melted_df)
pair <- data$interacting_pair

library(dplyr)
data <- read.table('/home/jhy/single_cell_TNBC/result/LR_241004/output/case/statistical_analysis_means_10_07_2024_110111.txt',sep='\t',header=T)
data <- subset(data, data$interacting_pair %in% pair)

#data <- read.table('/home/jhy/single_cell_TNBC/result/LR_241004/output/case/case_p0.05_LR_simple_meanfile.txt',sep='\t',header=T)
#cells <- colnames(data)[14:29]
cells <- colnames(data)[14:22]
cells <- lapply(cells, function(x) gsub('cell.','cell - ',x))
cells <- lapply(cells, function(x) gsub('\\.', ' ', x))
colnames(data) <- c(colnames(data)[1:13],cells)
create_combinations <- function(list_a) {
        result <- character()
        for (i in seq_along(list_a)) {
                for (j in seq_along(list_a)) {
                        result <- c(result, paste0(list_a[i], " - ", list_a[j]))
                }
        }
        return(result)
}
#a <- c('B cell','Myeloid cell','NK cell','T cell')
a <- c('Central memory T cell','Effector T cell','Regulatory T cell')
b <- create_combinations(a)
colname <- c(colnames(data)[1:13],b)
new <- data[,colname]
#test <- t(new[,14:29])
test <- t(new[,14:22])
colnames(test) <- new$interacting_pair
melted_df2 <- melt(test, id.vars = rownames(test))
colnames(melted_df2) <- c("cell_cell", "gene_gene", "mean")
head(melted_df2)
melted_df$mean <- melted_df2$mean
head(melted_df)

data <- melted_df

split_cell_cell <- strsplit(as.character(data$cell_cell), " - ")
data$cell1 <- sapply(split_cell_cell, `[`, 1)
data$cell2 <- sapply(split_cell_cell, `[`, 2)
#order_sequence <- c('T cell','NK cell','B cell','Myeloid cell')

#data$cell3 <- as.numeric(factor(data$cell1, levels = order_sequence))
#data$cell_cell <- paste(data$cell2, data$cell3, sep='   ')
data$cell_cell <- paste0(data$cell1, '-', data$cell2, sep='   ')

unique_values <- unique(data$cell_cell)
data$cell_cell <- factor(data$cell_cell, levels = unique_values)

neg_log_pval <- -log10(data$pval)
data$neglogp <- ifelse(is.infinite(neg_log_pval), 3, neg_log_pval)

### Visualization ###
library(ggplot2)
pdf('/home/jhy/single_cell_TNBC/result/LR_241004/output/case/dotplot_p0.05_threshold0.2_onlyLRandsimple_method2.pdf')
ggplot(data=data,aes(x=cell_cell, y=gene_gene, color=log2(mean), size=neglogp)) + geom_point() + scale_color_viridis_c(name = 'log2 (mean)') +cowplot::theme_cowplot() + theme(axis.line  = element_blank()) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab('') + theme(axis.ticks = element_blank())
dev.off()



