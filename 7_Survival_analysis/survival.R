#! /usr/bin/env Rscript
  
# Updated by Hye-Yeon Ju on Sep 10 2025.

##### Code Description #####
# This code performs survival analysis using TCGA public data.

##### Code Example #####
# Rscript
# /home/jhy/single_cell_TNBC/code/integration/survival/integration/survival/2ndRev_250902/Survival_FISH.R CCL5 ENSG00000271503 /home/jhy/single_cell_TNBC/result/survival_UseThis_250903/survival_250903_tnbc_IHC_FISH_mtestDef/

library(survival)
library(ggplot2)
library(survminer)
library(dplyr)
library(maxstat)
library(stringr)

args = commandArgs(trailingOnly = TRUE)
gene <- args[1]
id <- args[2]
outputdir <- args[3]

## Import RNAseq data
count <- read.table('/home/jhy/single_cell_TNBC/files/TCGA-BRCA.star_fpkm-uq.tsv',sep='\t',header=T,check.names=F)
# Check samples
ids <- colnames(count)[-1]
suffix <- sub(".*-", "", ids)
table(suffix)
sam01 <- ids[grepl("-01", ids)]
sam06 <- ids[grepl("-06", ids)]
sam01_clean <- sub("-[^-]+$", "", sam01)
sam06_clean <- sub("-[^-]+$", "", sam06)
sam06_clean
sam06_in_sam01_true <- sam06_clean[sam06_clean %in% sam01_clean]
sam06_in_sam01_true
# Select samples
ids_01 <- ids[grepl("-01[A-C]?$", ids)]
patient_id <- sub("-[^-]+$", "", ids_01)
priority <- c("01A", "01B", "01C")
selected <- sapply(unique(patient_id), function(pid) {
			   sam <- ids_01[patient_id == pid]
			   suffix <- sub(".*-", "", sam)
			   sam[which.min(match(suffix, priority))]
})
selected <- unname(selected)
count <- count[, c(1, which(colnames(count) %in% selected))]

## Import survival data
surv <- read.table('/home/jhy/single_cell_TNBC/files/TCGA-BRCA.survival.tsv',sep='\t',header=T)
rownames(count) <- count$Ensembl_ID
# Select samples
patient_ID <- sub("-[^-]+$", "",colnames(count))
count_with_surv <- count[,patient_ID %in% surv$X_PATIENT]

## Import clinical data (IHC) for filtering TNBC
clin <- read.table('/home/jhy/single_cell_TNBC/files/brca_tcga_clinical_data.tsv',sep='\t',header=T,quote='')
# IHC & FISH filtering
neg1 <- subset(clin, `ER.Status.By.IHC` == "Negative" & `PR.status.by.ihc` == "Negative")
length(unique(neg1$Patient.ID)) #220
neg2 <- subset(neg1, is.na(`HER2.fish.status`) | `HER2.fish.status` != "Positive")
length(unique(neg2$Patient.ID)) #210
ihc_status_keep <- c("Negative", "Equivocal", "Indeterminate", NA)
neg3 <- subset(
	       neg2,
	       (HER2.fish.status == "Negative" & IHC.HER2 %in% ihc_status_keep) |
		       (IHC.HER2 %in% c("Negative") & HER2.ihc.score %in% c("0","1",NA))
	       )

length(unique(neg3$Patient.ID)) #157
neg <- neg3

count_with_surv_tnbc <- count_with_surv[,patient_ID %in% neg$Patient.ID]

new_rownames <- sub("\\.[^.]*$", "", rownames(count_with_surv_tnbc))
input <- as.data.frame(t(count_with_surv_tnbc[new_rownames == id,]))

rownames(surv) <- surv$sample
input <- merge(surv, input, by='row.names')
input <- input[,c(1,3,5,6)]
colnames(input) <- c('id','status','time','exp')

## Import clinical data (age, stage) for multivariate analysis
clin <- read.table('/home/jhy/single_cell_TNBC/files/clinical.tsv',header=T,sep='\t',fill=T,stringsAsFactors=F,quote="")
patient_ID <- sub("-[^-]+$", "",input$id)
clin_filtered <- clin[clin$case_submitter_id %in% patient_ID, ]
#selected_cols <- c("case_submitter_id", "age_at_index", "ajcc_pathologic_stage", "ajcc_pathologic_t", "ajcc_pathologic_n", "ajcc_pathologic_m")
selected_cols <- c("case_submitter_id", "age_at_index", "ajcc_pathologic_stage")
clin_selected <- clin_filtered[, selected_cols]
clin_selected$age_at_index <- as.numeric(as.character(clin_selected$age_at_index))
clin_selected <- clin_selected %>% mutate(stage_raw = ajcc_pathologic_stage, stage_clean = toupper(trimws(gsub("'", "", stage_raw))))
clin_selected <- clin_selected %>% distinct(case_submitter_id, .keep_all = TRUE)
test <- clin_selected %>% mutate(stage_group = case_when(grepl("^STAGE\\s*IV", stage_clean)  ~ "Stage IV",grepl("^STAGE\\s*III", stage_clean) ~ "Stage III",grepl("^STAGE\\s*II", stage_clean)  ~ "Stage II",grepl("^STAGE\\s*I", stage_clean) ~ "Stage I",TRUE ~ NA_character_))
table(test$stage_group, useNA = "ifany")
clin_selected <- clin_selected %>% mutate(stage_group = case_when(grepl("^STAGE\\s*IV", stage_clean)  ~ "Stage IV",grepl("^STAGE\\s*III", stage_clean) ~ "Stage III",grepl("^STAGE\\s*II", stage_clean)  ~ "Stage II",grepl("^STAGE\\s*I", stage_clean) ~ "Stage I",TRUE ~ NA_character_)) %>% filter(!is.na(stage_group))

input$patient_id <- sub("-[^-]+$", "",input$id)
merge <- input %>% inner_join(clin_selected %>% select(case_submitter_id, stage_group, age_at_index),by = c("patient_id" = "case_submitter_id"))

merge <- merge %>% mutate(stage_group2 = case_when(stage_group %in% c("Stage I", "Stage II") ~ "Stage I-II",stage_group %in% c("Stage III", "Stage IV") ~ "Stage III-IV"))
merge$stage_group2 <- factor(merge$stage_group2, levels = c("Stage I-II", "Stage III-IV"))
merge$stage_group2 <- relevel(merge$stage_group2, ref = "Stage I-II")
input <- merge
colnames(input)[colnames(input) == "age_at_index"] <- 'Age'
colnames(input)[colnames(input) == "stage_group2"] <- 'Stage'

## Import ESTIMATE score for multivariate analysis
EST <- read.table('/home/jhy/single_cell_TNBC/files/breast_cancer_RNAseqV2.txt',header=T,sep='\t')
input <- input %>% mutate(id_short = sub(".$", "", id)) %>% left_join(EST %>% rename(id_short = ID), by = "id_short")
input$Tumor_purity <- cos(0.6049872018+0.0001467884*input$ESTIMATE_score)
head(input)

## Check covariates for multivariate analysis
cont_vars <- input[,c('Immune_score','Tumor_purity','Age')]
cor(cont_vars, method = "pearson")
cor(cont_vars, method = "spearman")
norm_test <- apply(cont_vars, 2, function(x) shapiro.test(x)$p.value)
norm_test

## Cutoff for survival analysis
#Maxstat
mstat=maxstat.test(Surv(time,status)~exp,data=input,smethod="LogRank",pmethod="condMC")
#mstat=maxstat.test(Surv(time,status)~exp,data=input,smethod="LogRank",pmethod="condMC",minprop=0.3,maxprop=0.7)
cutpoint=mstat$estimate
cutpoint
#Median
#cutpoint <- median(input$exp)
#Surv_cutpoint
#cutpoint <- surv_cutpoint(input,time = "time", event = "status", variables = "exp")
#cutpoint <- surv_cutpoint(input,time = "time", event = "status", variables = "exp",minprop=0.1)
#cutpoint <- surv_cutpoint(input,time = "time", event = "status", variables = "exp",minprop=0.2)
#cutpoint <- surv_cutpoint(input,time = "time", event = "status", variables = "exp",minprop=0.3)

## Survival analysis
input[[gene]]=ifelse(input$exp<=cutpoint,'low','high')
input[[gene]]<-factor(input[[gene]])
#fit <- survfit(as.formula(paste("Surv(time, status) ~", gene)), data = input)
#Only for HLA-B and MT-ND4L
fit <- survfit(as.formula(paste0("Surv(time, status) ~ `", gene, "`")), data = input)
fit

#Univariate KM
pdf(paste0(outputdir,gene,"_Kaplan_Meier.pdf"))
ggsurvplot(fit,legend.title=paste0(gene,"_Expression"),conf.int=T,pval=T,surv.median.line="hv",risk.table=T)
dev.off()
png(paste0(outputdir,gene,"_Kaplan_Meier.png"))
ggsurvplot(fit,legend.title=paste0(gene,"_Expression"),conf.int=T,pval=T,surv.median.line="hv",risk.table=T)
dev.off()

input[[gene]] = relevel(input[[gene]], ref='low')
vars <- c(gene, "Age", "Stage", "Tumor_purity", "Immune_score")

#Univariate coxph
pdf(file = paste0(outputdir,gene,"_Hazard_ratio_uni.pdf"), width = 6, height = 4)
for (v in vars) {
	formula <- as.formula(paste("Surv(time, status) ~", v))
	fit <- coxph(formula, data = input)
	print(fit)
	print(summary(fit))
	p <- ggforest(fit, data = input)
	print(p)
}
dev.off()

#Multivariate coxph with immune score
#fit.coxph <- coxph(as.formula(paste("Surv(time, status) ~", gene, "+ Age + Stage + Tumor_purity + Immune_score")),data = input)
#Only for HLA-B and MT-ND4L
fit.coxph <- coxph(as.formula(paste0("Surv(time, status) ~ `", gene, "` + Age + Stage + Tumor_purity + Immune_score")),data = input)
fit.coxph
print(summary(fit.coxph))
pdf(paste0(outputdir,gene,"_Hazard_ratio_multi_with_immune.pdf"),width=14,height=7)
ggforest(fit.coxph, data = input)
dev.off()

#Multivariate coxph without immune score
#fit.coxph <- coxph(as.formula(paste("Surv(time, status) ~", gene, "+ Age + Stage + Tumor_purity")),data = input)
#Only for HLA-B and MT-ND4L
fit.coxph <- coxph(as.formula(paste0("Surv(time, status) ~ `", gene, "` + Age + Stage + Tumor_purity")),data = input)
fit.coxph
print(summary(fit.coxph))
pdf(paste0(outputdir,gene,"_Hazard_ratio_multi.pdf"),width=14,height=7)
ggforest(fit.coxph, data = input)
dev.off()


