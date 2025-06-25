#!/bin/bash
#SBATCH -J gapit6 --mem 100G
#SBATCH -N 1 # number of nodes
#SBATCH -c 1 # number of tasks
#SBATCH -p ei-long
#SBATCH --mail-type=FAIL # notifications for job done & fail
#SBATCH --mail-user=kate.denning-james@earlham.ac.uk #send-to address

singularity exec /ei/projects/3/3c493fe7-33a3-4d92-b57f-50d4e29c3391/scratch/KDJ_beans/bwa_Andean/gatk_filter_GWAS/GAPIT_chrs/pca3_bytrait/gapit_r_4_2_1.v2.img Rscript --vanilla --verbose gapit_pca3_thin5_bytrait.R 
