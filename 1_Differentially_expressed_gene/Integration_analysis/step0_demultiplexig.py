#! /usr/bin/python2

# Updated by Hye-Yeon Ju on Jun 19, 2023.

##### Code Description #####
# This code performs demultiplexing of feature barcoded samples (TNBC and healthy donor data).
# You need to get the latest version of cellranger.
# Before run this code, you need to make csv files with reference path, fastq path, sample id, library type, and hash tags and save it to the result directory. Please check the format of csv files in resultdir.

##### Code Example #####
# python_script + raw_data_dir + resultdir
# /home/jhy/single_cell_TNBC/code/for_github/1_Differentially_expressed_gene/Integration_analysis/step0_demultiplexig.py /data1/analysis/BC_INU_VU_KU/03_TNBC_metastasis_scRNAseq/BCG_analysis/02_Spliting_samples/

import sys,os
import parameter as pr

resultdir = sys.argv[1]
sampledirs = ['TK-1','TK-2']

### Demultiplexing case and control samples ###
job="cellranger"
for sample in sampledirs:
    csv = resultdir + '/5p_hashing_demux_'+sample+'_config.csv'
    command = 'sbatch -J ' + job + ' -D ' + resultdir + ' --wrap="cellranger multi --id=demultiplexed_samples --csv=' + csv + '"'
    #print command
    pr.inShell(command)
    time.sleep(50)

