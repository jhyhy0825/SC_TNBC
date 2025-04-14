# SC_TNBC

##################################################
## Single-cell RNA sequencing analysis pipeline ##
##################################################

Updated by Hye-Yeon Ju on Apr 14, 2025.

- This pipeline is for single-cell RNA sequencing analysis using TNBC blood datasets.
- Each code includes 'Code Description' and 'Code Example'.
- Run 1_Differentially_expressed_gene > 2_Ligand_receptor_crosstalk > ... > 7_Survival_analysis
- Run step1_* > step2_* > ... > stepn_* scripts in each folder.

<<< 1_Differentially_expressed_gene >>>
Gene expression estimation, dimensionality reduction, clustering, DEG analysis, and visualization of scRNA-seq data using cellranger(10X genomics) and Seurat(R package).
Integration_analysis folder contains pipeline for integration analysis using internal and public datasets. Run step1 ~ step5 script.
Public_data_validation folder contains pipeline for validation using public datasets (gene expression profiling using pre vs post chemotherapy patients, chemotherapy responder vs non-responder patients)
Public datasets information: GSE169246 (mTNBC blood samples) / GSE174609, GSE134355 (healthy donor blood samples)

<<< 2_Ligand_receptor_crosstalk >>>
Ligand-receptor crosstalk analysis using results from 1_Differentially_expressed_gene.

<<< 3_Pseudotime_trajectory >>>
Pseudotime trajectory analysis using results from 1_Differentially_expressed_gene.

<<< 4_High_dimensional_WGCNA >>>
High dimensional weighted gene co-expression network analysis using results from 1_Differentially_expressed_gene.

<<< 5_Copy_number_variation >>>
Copy number variation analysis using results from 1_Differentially_expressed_gene.

<<< 6_Pathway_enrichment_analysis >>>
Pathway enrichment analysis using results from 1_Differentially_expressed_gene.
This code performs GO and KEGG pathway enrichment analysis.

<<< 7_Survival_analysis >>>
This code performs overall survival analysis using TCGA TNBC public dataset.
This code visualizes Kaplan-Meier plot and calculates hazard ratio of each gene.

