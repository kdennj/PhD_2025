#!/bin/bash
#SBATCH -p ei-short
#SBATCH -c 1 -N 1 --mem 20G
#SBATCH -J bedgraph_med
#SBATCH --mail-user=denning@nbi.ac.uk
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --array=0-144%15

cd /hpc-home/denning/KDJ_Scratch_project/KDJ_beans/bowtie2_meso_andean/unique/bedgraphs
REFDIR=/hpc-home/denning/KDJ_Scratch_project/KDJ_beans/ref_genomes/Andean_Meso_ref_genome

mapfile -t file < bedgraphs.txt

source bedtools-2.24.0

bedtools map -a ${REFDIR}/Andean_Meso_ref.fa.fastalength.txt.100kbwindows.txt -b ${file[${SLURM_ARRAY_TASK_ID}]} -c 4 -o median > median_cov/${file[${SLURM_ARRAY_TASK_ID}]}_median_cov.txt
