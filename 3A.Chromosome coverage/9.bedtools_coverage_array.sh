#!/bin/bash
#SBATCH -p ei-short
#SBATCH -c 1 -N 1 --mem 20G
#SBATCH -J bedtools_cov
#SBATCH --mail-user=denning@nbi.ac.uk
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --array=0-143%20

cd /KDJ_beans/bowtie2_meso_andean/unique

REFDIR=/KDJ_beans/ref_genomes/Andean_Meso_ref_genome

mapfile -t file < bam.txt

source bedtools-2.24.0

genomeCoverageBed -bg -split -ibam ${file[${SLURM_ARRAY_TASK_ID}]} -g ${REFDIR}/AM_window100kb_rep.bed > bedgraphs/${file[${SLURM_ARRAY_TASK_ID}]}.bedgraph