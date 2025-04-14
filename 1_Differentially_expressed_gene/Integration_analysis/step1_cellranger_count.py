#! /usr/bin/python2

# Updated by Hye-Yeon Ju on Jun 19, 2023.

##### Code Description #####
# This code performs mapping to the reference genome and gene expression estimation.
# You need to get the latest version of cellranger and reference data from 10X Genomics. 
# Before run this code, you need to make library files with fastq path, sample id, and lirary type.

##### Code Example #####
# python_script + raw_data_dir + resultdir
# /home/jhy/single_cell_TNBC/code/for_github/1_Differentially_expressed_gene/Integration_analysis/step1_cellranger_count.py /data1/analysis/BC_INU_VU_KU/03_TNBC_metastasis_scRNAseq/BCG_analysis/03_bamtofastq/ /data1/analysis/BC_INU_VU_KU/03_TNBC_metastasis_scRNAseq/BCG_analysis/04_gene_expression_profiling/

import sys,os
import parameter as pr

inputdir = sys.argv[1]
resultdir = sys.argv[2]
refdir = '/home/jhy/test/230217_scRNAseq/ref/refdata-gex-GRCh38-2020-A/'
feature = '/home/jhy/single_cell_TNBC/files/feature_reference_TotalSeqC_and_FB.csv'
library = '/home/jhy/single_cell_TNBC/files/'

### Gene expression counting ###
sampledirs = os.listdir(inputdir)
job="cellranger"
for sample in sampledirs:
    sample_id = sample+'-Control'
    sample_library = library+'libraries_'+sample+'-Control.csv'
    sample_result = resultdir+sample_id+'/'
    pr.makedir(sample_result)
    command1 = 'sbatch -J ' + job + ' -D ' + sample_result + ' --wrap="cellranger count --id ' + sample_id + ' --transcriptome=' + refdir + ' --libraries=' + sample_library + ' --feature-ref=' + feature + '"'
    #print command1
    pr.inShell(command1)
    sample_id = sample+'-Case'
    sample_library = library+'libraries_'+sample+'-Case.csv'
    sample_result = resultdir+sample_id+'/'
    pr.makedir(sample_result)
    command2 = 'sbatch -J ' + job + ' -D ' + sample_result + ' --wrap="cellranger count --id ' + sample_id + ' --transcriptome=' + refdir + ' --libraries=' + sample_library + ' --feature-ref=' + feature + '"'
    #print command2
    pr.inShell(command1)
    time.sleep(50)

