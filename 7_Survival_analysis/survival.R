#! /usr/bin/env Rscript
  
# Updated by Hye-Yeon Ju on Aug 7 2024.

##### Code Description #####
# This code performs survival analysis using TCGA public data.

##### Code Example #####
# Rscript
# /home/jhy/single_cell_TNBC/code/integration/survival/survival.R

id <- 'ENSG00000213402'
gene <- 'PTPRCAP'

### Import input count data ###
count1 <- read.table('/home/jhy/single_cell_TNBC/files/TCGA-BRCA.htseq_counts.tsv',sep='\t',header=T)
#count1 <- read.table('/home/jhy/single_cell_TNBC/files/TCGA-BRCA.htseq_fpkm-uq.tsv',sep='\t',header=T)
count2 <- count1[, !grepl("\\.11", colnames(count1))]
###count2 <- count2[, !grepl("\\.06", colnames(count2))]
#last_parts <- sapply(strsplit(colnames(count2), "\\."), function(x) tail(x, 1))
#unique_last_parts <- unique(last_parts)
#unique_last_parts
#length(colnames(count2))
#test1 <- colnames(count2)

### Import input survival data ###
surv1 <- read.table('/home/jhy/single_cell_TNBC/files/TCGA-BRCA.survival.tsv',sep='\t',header=T)
#print(unique_values)
surv2 <- surv1[!grepl("-11", surv1$sample), ]
###surv2 <- surv2[!grepl("-06", surv2$sample), ]
#extracted_values <- sub(".*-", "", surv2$sample)
#unique_values <- unique(extracted_values)
#unique_values
#dim(surv2)
test2 <- gsub("-", ".", surv2$sample)

#not_in_test2 <- test1[!test1 %in% test2]
#not_in_test1 <- test2[!test2 %in% test1]

#gdc <- read.table('/home/jhy/single_cell_TNBC/files/clinical.tsv',sep='\t',header=T,fill=TRUE,quote="")
#ids <- gdc$case_submitter_id
#test2a <- sub("\\.[^.]*$", "", test2)
#test3 <- gsub("-", ".", ids)
#not_in_test2a <- test3[!test3 %in% test2a]
#not_in_test3 <- test2a[!test2a %in% test3]
#gdc$new_id <- test3
#not_in_test2a <- sub("\\.[^.]*$", "", not_in_test2)
#new <- subset(gdc, gdc$new_id %in% not_in_test2a)
#new$days_to_death

rownames(count2) <- count2$Ensembl_ID
new <- count2[,colnames(count2) %in% test2]
new_rownames <- sub("\\.[^.]*$", "", rownames(new))

#PTPRCAP <- as.data.frame(t(new[new_rownames == 'ENSG00000213402',]))
#HOPX <- as.data.frame(t(new[new_rownames == 'ENSG00000171476',]))
#GZMK <- as.data.frame(t(new[new_rownames == 'ENSG00000113088',]))
#CTSW <- as.data.frame(t(new[new_rownames == 'ENSG00000172543',]))
#KLRG1 <- as.data.frame(t(new[new_rownames == 'ENSG00000139187',]))
#PRF1 <- as.data.frame(t(new[new_rownames == 'ENSG00000180644',]))
#NKG7 <- as.data.frame(t(new[new_rownames == 'ENSG00000105374',]))
#JUNB <- as.data.frame(t(new[new_rownames == 'ENSG00000171223',]))
#CCL4 <- as.data.frame(t(new[new_rownames == 'ENSG00000275302',]))
#GNLY <- as.data.frame(t(new[new_rownames == 'ENSG00000115523',]))
pam <- read.table('/home/jhy/single_cell_TNBC/files/BRCA.547.PAM50.SigClust.Subtypes.txt',header=T,sep='\t')
pam <- subset(pam, pam$PAM50 == 'Basal')
pam$Sample <- sapply(pam$Sample, function(x) paste(unlist(strsplit(x, "-"))[1:4], collapse = "-"))
pam$Sample <- gsub("-", ".", pam$Sample)
outputdir <- '/home/jhy/single_cell_TNBC/result/Seurat/integration2/result/survival/TNBC_count_max/'
library(survival)
library(ggplot2)
library(survminer)
library(dplyr)
library(maxstat)

input <- as.data.frame(t(new[new_rownames == id,]))
input <- subset(input,rownames(input) %in% pam$Sample)

surv2$sample <- gsub("-", ".", surv2$sample)
rownames(surv2) <- surv2$sample

input <- merge(surv2, input, by='row.names')
input <- input[,c(1,3,5,6)]
colnames(input) <- c('id','status','time','exp')

mstat=maxstat.test(Surv(time,status)~exp,data=input,smethod="LogRank",pmethod="condMC")
cutpoint=mstat$estimate
cutpoint

#cutpoint <- mean(input$exp)
input$group=ifelse(input$exp<=cutpoint,'low','high')
#input$group=ifelse(input$exp<cutpoint,'low','high')
input$group<-factor(input$group)
fit<-survfit(Surv(time,status)~group, data=input)
fit

### Visualization ###
pdf(paste0(outputdir,gene,"_Kaplan_Meier.pdf"))
ggsurvplot(fit,legend.title=paste0(gene,"_Expression"),conf.int=T,pval=T,surv.median.line="hv",risk.table=T)
dev.off()

png(paste0(outputdir,gene,"_Kaplan_Meier.png"))
ggsurvplot(fit,legend.title=paste0(gene,"_Expression"),conf.int=T,pval=T,surv.median.line="hv",risk.table=T)
dev.off()

input$group = relevel(input$group, ref='low')
fit.coxph <- coxph(Surv(time,status)~group, data = input)
fit.coxph
pdf(paste0(outputdir,gene,"_Hazard_ratio.pdf"),width=14,height=7)
ggforest(fit.coxph, data = input)
dev.off()


