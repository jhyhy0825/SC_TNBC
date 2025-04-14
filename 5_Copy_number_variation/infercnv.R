#! /usr/bin/env Rscript
  
# Updated by Hye-Yeon Ju on May 24, 2024.

##### Code Description #####
# This code performs copy number variation analysis using infercnv.

##### Code Example #####
# Rscript
# /home/jhy/single_cell_TNBC/code/for_github/5_Copy_number_variation/infercnv.R

library(Seurat)
library(infercnv)

### Import RDS file ###
data <- readRDS('/home/jhy/single_cell_TNBC/result/Seurat/integration2/result2/Integration_mt10_feat4000_dim10_res0.4_with_case_control.rds')
DefaultAssay(object = data) <- "SCT"
counts_matrix = data[["SCT"]]$counts
anno <- data$orig.ident
write.table(counts_matrix,'/home/jhy/single_cell_TNBC/result/Seurat/integration2/result2/test_count.tsv',sep='\t', quote=F)
list2 <- names(anno)

### It depends on the sample type ###
#list2_modified <- lapply(list2, function(x) {
# modified_string <- gsub("-", ".", x)
#  return(paste0("X", modified_string))
#})

###
list2_modified <- lapply(list2, function(x) {
  modified_string <- gsub("-", ".", x)
  if (grepl("^8674", modified_string)) {
    return(paste0("X", modified_string))
  } else {
    return(modified_string)
  }
})
###

### Run infercnv ###
names(anno) <- list2_modified
write.table(anno,'/home/jhy/single_cell_TNBC/result/Seurat/integration2/result2/test_annotation.tsv',sep='\t',quote=F,col.names=F)
infercnv_obj = CreateInfercnvObject(raw_counts_matrix='/home/jhy/single_cell_TNBC/result/Seurat/integration2/result2/test_count.tsv',gene_order_file='/home/jhy/single_cell_TNBC/files/hg38_gencode_v27.txt',annotations_file='/home/jhy/single_cell_TNBC/result/Seurat/integration2/result2/test_annotation.tsv', delim='\t', ref_group_names=NULL)
infercnv_obj = infercnv::run(infercnv_obj, cutoff=1, out_dir='/home/jhy/single_cell_TNBC/result/Seurat/integration2/result2/infercnv/', cluster_by_groups=TRUE, denoise=TRUE, HMM=TRUE,output_format='pdf')

### Visualization ###
plot_cnv(infercnv_obj,out_dir = '/home/jhy/single_cell_TNBC/result/Seurat/integration2/result2/infercnv4/',write_expr_matrix = TRUE, cluster_by_groups=TRUE,output_format='pdf')

mean(as.numeric(unlist(data['KLRD1',])),na.rm = TRUE)
