#! /usr/bin/env Rscript
  
# Updated by Hye-Yeon Ju on Apr 23, 2024.

##### Code Description #####
# This code performs pathway enrichment analysis using clusterProfiler.

##### Code Example #####
# Rscript
# /home/jhy/single_cell_TNBC/code/for_github/6_Pathway_enrichment_analysis/clusterProfiler.R

library(dplyr)
library(Seurat)
library(patchwork)
library("clusterProfiler")
library("org.Hs.eg.db")
library("AnnotationHub")
library(DOSE)
library(enrichplot)
library(dplyr)
library(biomaRt)

### Import RDS file ###
data <- readRDS('/home/jhy/single_cell_TNBC/result/Seurat/integration2/result2/nkcell.rds')
data <- PrepSCTFindMarkers(data)
markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top100 <- markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
top100pval <- subset(top100, rowSums(top100[5] < 0.05) > 0)
test <- subset(top100pval, cluster=='0')
test <- data.frame(test)
adjP <- test$p_val_adj
gene <- test$gene
gene <- bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
gene <- gene$ENTREZID
de <- gene

### GO analysis ###
ego <- enrichGO(gene=de,OrgDb=org.Hs.eg.db)
pdf('/home/jhy/single_cell_TNBC/result/Seurat/integration2/result2/nkcell_barplot_go_cluster0.pdf',height=8,width=6)
barplot(ego, showCategory=20)
dev.off()

### KEGG analysis ###
ekg <- enrichKEGG(gene=gene,organism='hsa')
pdf('/home/jhy/single_cell_TNBC/result/Seurat/integration2/result2/nkcell_barplot_kegg_cluster0.pdf',height=8,width=6)
barplot(ekg, showCategory=20)
dev.off()


## UseThis and correct cell type ##
clu='7'

test <- subset(top100pval, cluster==clu)
test <- data.frame(test)
adjP <- test$p_val_adj
gene <- test$gene
gene <- bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
gene <- gene$ENTREZID
de <- gene

ego <- enrichGO(gene=de,OrgDb=org.Hs.eg.db)
pdf(paste0('/home/jhy/single_cell_TNBC/result/Seurat/integration2/result2/bcell_barplot_go_cluster',clu,'.pdf'),height=8,width=6)
barplot(ego, showCategory=20)
dev.off()

ekg <- enrichKEGG(gene=gene,organism='hsa')
pdf(paste0('/home/jhy/single_cell_TNBC/result/Seurat/integration2/result2/bcell_barplot_kegg_cluster',clu,'.pdf'),height=8,width=6)
barplot(ekg, showCategory=20)
dev.off()


