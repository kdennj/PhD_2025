#!/bin/bash
#SBATCH -p ei-short
#SBATCH -c 4 -N 1 --mem 20G
#SBATCH -J bedgraph_mean
#SBATCH --mail-user=denning@nbi.ac.uk
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --array=0-144%15

cd /KDJ_beans/bowtie2_meso_andean/unique/bedgraphs
REFDIR=/KDJ_beans/ref_genomes/Andean_Meso_ref_genome

mapfile -t file < bedgraphs.txt

source bedtools-2.24.0

bedtools map -a ${REFDIR}/Andean_Meso_ref.fa.fastalength.txt.100kbwindows.txt \
-b  ${file[${SLURM_ARRAY_TASK_ID}]} -c 4 -o mean > mean_cov/${file[${SLURM_ARRAY_TASK_ID}]}_mean_cov.txt
